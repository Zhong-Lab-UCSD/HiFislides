#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    ---------------------------------------------------------------------------
    Usage: $PROGNAME [-b <bin_dir>] [-i <star_index>] [-g <ref_gtf>]
                     [-N <sample_name>]
                     [-S <spatial_barcode_dir>] 
                     [-1 <fastq.gz_R1>] [-2 <fastq.gz_R2>] 
                     [-t <threads>] [-o <output_dir>] 
    
    Dependency: pear, fastp, bwa, star, bedtools, samtools, seqtk
    This is an all-in-one wrapper script of HiFi Pipeline.
    -b : Directory of the scripts.
    -i : Directory of the STAR index.
    -g : GTF annotation file of the reference genome.
    -N : Sample name, used to label the final output files.
    -S : Directory of the processed spatial barcodes.
    -1 : R1 fastq.gz file of the HiFi library.
    -2 : R2 fastq.gz file of the HiFi library.
    -t : Max CPU threads for parallelized processing, at least 4 (default 8).
    -o : Parent output directory.
    -h : Show usage help
    ---------------------------------------------------------------------------
EOF
    exit 1
}

parameter_error() {
    echo "!!!!!!!!!!!!!!!! Pipeline exit with parameter ERROR at $(date '+%Y-%m-%d %H:%M:%S %Z') !!!!!!!!!!!!!!!!"
    usage
}

while getopts :b:i:g:N:S:1:2:t:o:h opt; do
    case $opt in
        b) BIN_DIR=${OPTARG};;
        i) STAR_INDEX=${OPTARG};;
        g) annotation_gtf_file=${OPTARG};;
        N) SAMPLE_NAME=${OPTARG};;
        S) L1_DIR=${OPTARG};;
        1) L2R1_FASTQ=${OPTARG};;
        2) L2R2_FASTQ=${OPTARG};;
        t) N_THREADS=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        h) usage;;
    esac
done

# check parameters
[ -z "$BIN_DIR" ] &&  echo "Error!! Please provide path to scripts with -b" && parameter_error

[ -z "$STAR_INDEX" ] && echo "Error!! Please provide path to the STAR index files with -i" && parameter_error

[  -z "$annotation_gtf_file" ] && echo "Error!! Please provide the GTF annotation file with -g" && parameter_error

[  -z "$SAMPLE_NAME" ] && echo "Error!! Please provide the sample name with -N" && parameter_error

[ ! -d "$L1_DIR" ] && echo "Error!! Directory of spatial barcodes does not exist: "$L1_DIR && parameter_error

[ -z "$L2R1_FASTQ" ] && echo "Error!! Please provide fastq.gz files of R1 with -1" && parameter_error
[ -z "$L2R2_FASTQ" ] && echo "Error!! Please provide fastq.gz files of R2 with -2" && parameter_error

[  -z "$N_THREADS" ] && echo "Use default thread number 8'." && N_THREADS=8
if ! [[ "$N_THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && parameter_error 
fi

[ ! -d "$OUT_DIR" ] && echo "Error!! Output directory not exist: "$OUT_DIR && parameter_error


################## INITIALIZE VARIABLES

L2_DIR=$OUT_DIR/$SAMPLE_NAME # HiFi library
mkdir -p $L2_DIR

L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping
mkdir -p $L2R2_GENOME_DIR

L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping
mkdir -p $L2R1_MAPPING_DIR

# L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate_N1
# mkdir -p $L2R1_L2R2_INTEGRATE_DIR

# Spatial barcode files
FLOWCELL_FULL=$(basename $L1_DIR)
L1R1_DUP=$L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt # second output of surfdedup
SEQ_MACHINE_ID=$(cut -f1 $L1R1_DUP | head -1 | cut -f1 -d ":")

L1R1_FASTQ_BWA_INDEX=$L1_DIR/bwa_index_L1R1/$FLOWCELL_FULL.L1R1_dedup # bwa index path and basename
L1R1_FASTQ_BWA_INDEX_TOP=$L1_DIR/bwa_index_L1R1_top/$FLOWCELL_FULL.L1R1_dedup.top # bwa index path and basename
L1R1_FASTQ_BWA_INDEX_BOTTOM=$L1_DIR/bwa_index_L1R1_bottom/$FLOWCELL_FULL.L1R1_dedup.bottom # bwa index path and basename

# Select full gene coordinates only
annotation_gtf_file_genes="$(dirname "${annotation_gtf_file}")"/gencode.v41.annotation.gene.gtf
awk -v OFS='\t' '$3=="gene"' $annotation_gtf_file > $annotation_gtf_file_genes


################## PROCESSING
> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

START_DATE=$(date) # start processing date
echo "Processing of "$SAMPLE_NAME
echo "Processing of "$SAMPLE_NAME >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


#################### LIBRARY 2 R2
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L2R2..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L2R2..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Preprocessing of HiFi R2 reads
mkdir -p $L2_DIR/L2R2_preprocessing

# Find L2R2 reads overlapping L2R1 and filter them out using the software PEAR
echo "[$(date '+%m-%d-%y %H:%M:%S')] Find L2R2 reads overlapping L2R1 and filter them out using the software PEAR..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

minoverlap=10 # default for PEAR

pear \
-f $L2R1_FASTQ \
-r $L2R2_FASTQ \
-v $minoverlap \
-j $N_THREADS \
-o $L2_DIR/L2R2_preprocessing/L2R2_pear

size=$(stat -c %s $L2_DIR"/L2R2_preprocessing/L2R2_pear.unassembled.reverse.fastq")

if [ $size == 0 ]; then
echo "[$(date '+%m-%d-%y %H:%M:%S')] No output file from PEAR. Stopping!"
echo "[$(date '+%m-%d-%y %H:%M:%S')] No output file from PEAR. Stopping!" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
exit 1
fi

seqtk seq -r $L2_DIR/L2R2_preprocessing/L2R2_pear.unassembled.reverse.fastq > $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq
echo "[$(date '+%m-%d-%y %H:%M:%S')] PEAR processing complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

rm $L2_DIR/L2R2_preprocessing/L2R2_pear* # FILES ARE BIG!!!

echo "[$(date '+%m-%d-%y %H:%M:%S')] Filtering L2R2 reads..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

if [ $N_THREADS -ge 16 ]; then
fastp_thread=16
else
fastp_thread=$N_THREADS
fi

### FASTP filtering
fastp \
-i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
-o $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.fastq \
-h $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.log.html \
-j $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.log.json \
--trim_poly_g \
--trim_poly_x \
--disable_quality_filtering \
--thread $fastp_thread

echo "[$(date '+%m-%d-%y %H:%M:%S')] Filtering complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Align HiFi R2 reads to genome/genes in order to obtain gene annotation for HiFi read pairs.
echo "[$(date '+%m-%d-%y %H:%M:%S')] Align HiFi-Slide reads R2 (L2R2) to the genome using STAR..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

STAR \
--genomeDir $STAR_INDEX \
--readFilesIn $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.fastq \
--outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx \
--outSAMattributes All \
--outFileNamePrefix $L2R2_GENOME_DIR/L2R2_genome. \
--sjdbGTFfile $annotation_gtf_file \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--runThreadN $N_THREADS

### Select uniquely mapped reads
samtools view -@ $N_THREADS -b -h -q 255 \
-o $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.out.bam \
$L2R2_GENOME_DIR/L2R2_genome.Aligned.out.bam

### Map uniquely mapped reads over genes. Cannot use featureCounts because we need to keep track of what L2R2 read align to each gene.
bedtools intersect \
-a $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.out.bam \
-b $annotation_gtf_file_genes \
-wb -bed | cut -f 4,21 > $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp.bed

cut -f 2 $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp.bed |  # extract column 2
cut -f1 -d';' |        # extract column 1 using semi-colon as delimiter
cut -f2 -d' ' |        # extract column 2 using space as delimiter
tr -d '"' > $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_id.txt  # remove double quote

cut -f 2 $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp.bed |  # extract column 2
cut -f3 -d';' |        # extract column 3 using semi-colon as delimiter
cut -f3 -d' ' |        # extract column 3 using space as delimiter
tr -d '"' > $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_name.txt  # remove double quote

cut -f 2 $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp.bed |  # extract column 2
cut -f2 -d';' |        # extract column 2 (column 5 for Ensembl) using semi-colon as delimiter
cut -f3 -d' ' |        # extract column 3 using space as delimiter
tr -d '"' > $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_biotype.txt  # remove double quote

paste \
$L2R2_GENOME_DIR/HiFi_L2R2_genome_temp.bed \
$L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_id.txt \
$L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_name.txt \
$L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_biotype.txt | cut -f 1,3,4,5 > $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_1.bed

rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp.bed
rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_id.txt
rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_name.txt
rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_biotype.txt

### Reads not overlapping annotated genes 
bedtools intersect \
-a $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.out.bam \
-b $annotation_gtf_file_genes \
-v -bed | cut -f 4 | awk -v OFS='\t' '{print $1, "NA", "NA", "NA"}' > $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed

### Concatenate and sort
cat $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_1.bed $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed

rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_1.bed
rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed

# sed 's/ //g' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed > $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL1.sort.bed

echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment of HiFi-Slide reads R2 (L2R2) to the genome complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide library 2 finished."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide library 2 finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


#################### LIBRARY 2 R1
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L2R1..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L2R1..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


###### Align HiFi R1 reads (L2R1) to spatial barcodes (L1R1) in order to obtain spatial coordinates for HiFi read pairs.
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start aligning HiFi-Slide R1 reads (L2R1) to spatial barcodes (L1R1)..." 
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start aligning HiFi-Slide R1 reads (L2R1) to spatial barcodes (L1R1)..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Alignment to all the barcodes (-k 19 is default in bwa)
bwa mem \
-a \
-k 19 \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX $L2R1_FASTQ 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.log | grep -P "\t0\t$SEQ_MACHINE_ID|\t256\t$SEQ_MACHINE_ID" | awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed - > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.filter.sam

### Alignment to the barcodes in the top-half region of the flowcell (-k 19 is default in bwa)
bwa mem \
-a \
-k 19 \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX_TOP $L2R1_FASTQ 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.top.log | grep -P "\t0\t$SEQ_MACHINE_ID|\t256\t$SEQ_MACHINE_ID" | awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed - > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.filter.sam

### Alignment to the barcodes in the bottom-half region of the flowcell (-k 19 is default in bwa)
bwa mem \
-a \
-k 19 \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX_BOTTOM $L2R1_FASTQ 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.bottom.log | grep -P "\t0\t$SEQ_MACHINE_ID|\t256\t$SEQ_MACHINE_ID" | awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed - > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottom.filter.sam


bwa mem \
-a \
-k 19 \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX_TOP $L2R1_FASTQ > $L2R1_MAPPING_DIR/top.sam 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.top.log

bwa mem \
-a \
-k 19 \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX_BOTTOM $L2R1_FASTQ > $L2R1_MAPPING_DIR/bottom.sam 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.bottom.log

grep -v '^@' $L2R1_MAPPING_DIR/top.sam | grep -P"\t0\t$SEQ_MACHINE_ID|\t256\t$SEQ_MACHINE_ID" > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.sam
grep -v '^@' $L2R1_MAPPING_DIR/bottom.sam | grep -P "\t0\t$SEQ_MACHINE_ID|\t256\t$SEQ_MACHINE_ID" > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottom.sam

grep -v '^@' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.toph.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.sam
grep -v '^@' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottomh.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottom.sam


samtools view -F 4095 -o sam0.sam $L2R1_MAPPING_DIR/bottom.sam 
samtools view -f 256 -o sam256.sam $L2R1_MAPPING_DIR/bottom.sam 
samtools view -F 4095 -f 256 -o sam0256.sam $L2R1_MAPPING_DIR/bottom.sam

samtools view -F 4095 $L2R1_MAPPING_DIR/bottom.sam | wc -l
samtools view -f 256 $L2R1_MAPPING_DIR/bottom.sam | wc -l
samtools view -F 4095 -f 256 $L2R1_MAPPING_DIR/bottom.sam | wc -l

echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment done."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment done." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

## To be tested for Cortex_4
bwa mem \
-a \
-k 19 \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX_TOP $L2R1_FASTQ 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.top.log | grep -P "\t0\t$SEQ_MACHINE_ID|\t256\t$SEQ_MACHINE_ID" | awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed - > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.filter.sam


grep -P "\t0\t$SEQ_MACHINE_ID|\t256\t$SEQ_MACHINE_ID" $L2R1_MAPPING_DIR/top.sam | awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed - > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.filter.sam

grep -P "\t0\t$SEQ_MACHINE_ID|\t256\t$SEQ_MACHINE_ID" $L2R1_MAPPING_DIR/bottom.sam | awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed - > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottom.filter.sam




### Filter SAM file to select only HiFi-Slide reads mapped to genome/transcriptome (samf: "SAM filter" custom format)
echo "[$(date '+%m-%d-%y %H:%M:%S')] Filter SAM file to select only HiFi-Slide reads mapped to genome/transcriptome..."

# awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.filter.sam

# awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.filter.sam

# awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottom.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottom.filter.sam

# rm $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.sam
# rm $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.sam
# rm $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottom.sam

echo "[$(date '+%m-%d-%y %H:%M:%S')] Filter SAM file to select only HiFi-Slide reads mapped to genome/transcriptome complete."

### Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location
echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location..."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

python3 $BIN_DIR/hifiwrangling0.py \
$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.filter.sam \
$L1R1_DUP \
1000 | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifiwrangling0.sort.o

python3 $BIN_DIR/hifiwrangling0.py \
$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.filter.sam \
$L1R1_DUP \
1 | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifiwrangling0.N1.sort.o

python3 $BIN_DIR/hifiwrangling0.py \
$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.filter.sam \
$L1R1_DUP \
1000 | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifiwrangling0.top.sort.o

python3 $BIN_DIR/hifiwrangling0.py \
$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottom.filter.sam \
$L1R1_DUP \
1000 | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifiwrangling0.bottom.sort.o

python3 $BIN_DIR/hifiwrangling0.py \
$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.top.filter.sam \
$L1R1_DUP \
1 | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifiwrangling0.N1.top.sort.o

python3 $BIN_DIR/hifiwrangling0.py \
$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.bottom.filter.sam \
$L1R1_DUP \
1 | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifiwrangling0.N1.bottom.sort.o

echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location complete."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide L2R1 complete."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide L2R1 complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


########## Integrate spatial coordinates and gene expression information
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

for i in "." ".N1." ".top." ".bottom." ".N1.top." ".N1.bottom."; do

L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate$i
mkdir -p $L2R1_L2R2_INTEGRATE_DIR

### Genome
join -1 1 -2 1 -t $'\t' $L2R1_MAPPING_DIR"/L2R1_L1R1.hifiwrangling0"$i"sort.o" $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed | cut -f 2,3,4,5,7,8,9 > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt > $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt

# Calculate the number of read pairs in each tile
awk -v OFS='\t' '{print $1}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt | awk '{A[$1]++}END{for(i in A)print i,A[i]}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_read_number_table.txt

# Calculate the number of spots in each tile
awk -v OFS='\t' '{print $1, $1"_"$2"_"$3}' $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk '!seen[$2]++' | awk '{A[$1]++}END{for(i in A)print i,A[i]}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_spot_number_table.txt

# TO BE TESTED AND COPIED TO THE MAIN WRAPPER
# # Calculate the number of genes in each tile
# grep -v -P "\tNA\tNA\tNA\t" $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1, $1"_"$5}' | awk '!seen[$2]++' | awk '{A[$1]++}END{for(i in A)print i,A[i]}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_gene_number_table.txt

# # Calculate the total gene expression in each tile
# grep -v -P "\tNA\tNA\tNA\t" $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1, $7}' | awk -F"\t" '{array[$1]+=$2} END { for (i in array) {print i"\t" array[i]}}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_gene_expression_table.txt

done

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

#####################

echo "Data processing started: "$START_DATE
echo "Data processing ended: "$(date)

echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo "Data processing started: "$START_DATE >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo "Data processing ended: "$(date) >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log



