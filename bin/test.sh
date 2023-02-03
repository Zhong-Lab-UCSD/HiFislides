#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    ---------------------------------------------------------------------------
    Usage: $PROGNAME [-b <bin_dir>] [-i <star_index>] [-g <ref_gtf>]
                     [-N <sample_name>]
                     [-F <flowcell_id>] [-f <flowcell_type>] 
                     [-S <spatial_barcode_dir>] 
                     [-1 <fastq.gz_R1>] [-2 <fastq.gz_R2>] 
                     [-t <threads>] [-o <output_dir>] 
    
    Dependency: pear, fastp, bwa, star, bedtools, samtools, seqtk
    This is an all-in-one wrapper script of HiFi Pipeline.
    -b : Directory of the scripts.
    -i : Directory of the STAR index.
    -g : GTF annotation file of the reference genome.
    -N : Sample name, used to label the final output files.
    -F : Flowcell ID.
    -f : Flowcell type (NextSeq, MiniSeq).
    -S : Parent directory of the processed spatial barcodes.
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

while getopts :b:i:g:N:F:f:S:1:2:t:o:h opt; do
    case $opt in
        b) BIN_DIR=${OPTARG};;
        i) STAR_INDEX=${OPTARG};;
        g) annotation_gtf_file=${OPTARG};;
        N) SAMPLE_NAME=${OPTARG};;
        F) flowcell=${OPTARG};;
        f) flowcell_type=${OPTARG};;
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

[  -z "$SAMPLE_NAME" ] && && echo "Error!! Please provide the sample name with -N" && parameter_error

[  -z "$flowcell" ] && && echo "Error!! Please provide the flowcell ID with -F" && parameter_error

[  -z "$flowcell_type" ] && && echo "Error!! Please provide the flowcell type (NextSeq, MiniSeq) with -f" && parameter_error

[ ! -d "$L1_DIR" ] && echo "Error!! Directory of spatial barcodes does not exist: "$L1_DIR && parameter_error

[ -z "$L2R1_FASTQ" ] && echo "Error!! Please provide fastq.gz files of R1 with -1" && parameter_error
[ -z "$L2R2_FASTQ" ] && echo "Error!! Please provide fastq.gz files of R2 with -2" && parameter_error

[  -z "$N_THREADS" ] && echo "Use default thread number 8'." && N_THREADS=8
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && parameter_error 
fi

[ ! -d "$OUT_DIR" ] && echo "Error!! Output directory not exist: "$OUT_DIR && parameter_error


##############################

L2_DIR=$OUT_DIR/$SAMPLE_NAME # HiFi library
mkdir -p $L2_DIR

if [ "$flowcell_type" == "NextSeq" ]; then
surface=$flowcell:1:1
elif [ "$flowcell_type" == "MiniSeq" ]; then
surface=$flowcell:1:
fi

# Spatial barcode files
L1R1_FASTQ_BWA_INDEX=$L1_DIR/$flowcell/bwa_index_L1R1/$flowcell.L1R1_dedup # bwa index path and basename
L1R1_DEDUP=$L1_DIR/$flowcell/$flowcell.L1R1_dedup.fasta # first output of surfdedup
L1R1_DUP=$L1_DIR/$flowcell/$flowcell.L1R1_dup.txt # second output of surfdedup

# Select full gene coordinates only
annotation_gtf_file_genes="$(dirname "${annotation_gtf_file}")"/gencode.v41.annotation.gene.gtf
awk -v OFS='\t' '$3=="gene"' $annotation_gtf_file > $annotation_gtf_file_genes


################## PROCESSING
touch $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

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

L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping
L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.fastq

mkdir -p $L2R2_GENOME_DIR
STAR \
--genomeDir $STAR_INDEX \
--readFilesIn $L2R2_FILTER_FASTQ \
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
-v -bed | cut -f 4 | awk '{print $0, "\tNA", "\tNA", "\tNA"}' > $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed

### Concatenate and sort
cat $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_1.bed $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed

rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_1.bed
rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed

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

### Alignment for all the barcodes
L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping
mkdir -p $L2R1_MAPPING_DIR
bwa_seed_length=19 # default

bwa mem \
-a \
-k $bwa_seed_length \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX $L2R1_FASTQ 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup_k$bwa_seed_length.log | grep -v '^@' | awk -F"\t" '$2 == "0" || $2 == "256" { print $0 }' > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup_k$bwa_seed_length.sam

echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment done."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment done." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


### Select flowcell surface if flowcell type is MiniSeq
if [ "$flowcell_type" == "MiniSeq" ]; then
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selecting the correct flowcell surface..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# Select reads coming from surface 1 and surface 2
grep $surface"1" $L2R1_MAPPING_DIR/L2R1_L1R1_dedup_k$bwa_seed_length.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface1.sam
grep $surface"2" $L2R1_MAPPING_DIR/L2R1_L1R1_dedup_k$bwa_seed_length.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface2.sam

# Choose which surface the tissue is based on the number of mapped reads
n_reads_surface_1=$(wc -l $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface1.sam | cut -d " " -f 1)
n_reads_surface_2=$(wc -l $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface2.sam | cut -d " " -f 1)

if [ $n_reads_surface_1 > $n_reads_surface_2 ]; then
L2R1_L1R1_SAM=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface1.sam
rm $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface2.sam
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selection of the surface done. Selected surface 1." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
else
L2R1_L1R1_SAM=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface2.sam
rm $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface1.sam
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selection of the surface done. Selected surface 2." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
fi

elif [ "$flowcell_type" == "NextSeq" ]; then

L2R1_L1R1_SAM=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup_k$bwa_seed_length.sam

fi


### Filter SAM file to select only HiFi-Slide reads mapped to genome/transcriptome (samf: "SAM filter" custom format)
L2R1_L1R1_SAM_FILTER=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.filter.sam

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed $L2R1_L1R1_SAM > $L2R1_L1R1_SAM_FILTER

# rm $L2R1_L1R1_SAM

### Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location
echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location..."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

python3 $BIN_DIR/hifiwrangling0.py \
$L2R1_L1R1_SAM_FILTER \
$L1_DIR/$flowcell/$flowcell.L1R1_dup.txt \
1000 | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifiwrangling0.sort.o

echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location complete."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score and match HiFi-Slide read pairs with spatial location complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide L2R1 complete."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide L2R1 complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


########## Integrate spatial coordinates and gene expression information
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate
mkdir -p $L2R1_L2R2_INTEGRATE_DIR

### Genome
HiFi_L2R2_genome=$L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed

join -1 1 -2 1 -t $'\t' $L2R1_MAPPING_DIR/L2R1_L1R1.hifiwrangling0.sort.o $HiFi_L2R2_genome | cut -f 2,3,4,5,6,7,8 > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt > $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt

# rm $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt

# Calculate the number of spots in each tile
awk -v OFS='\t' '{print $1, $1"_"$2"_"$3}' $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk '!seen[$2]++' | awk '{A[$1]++}END{for(i in A)print i,A[i]}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_spot_number_table.txt

# Calculate the number of genes in each tile
awk -v OFS='\t' '{print $1, $1"_"$5}' $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk '!seen[$2]++' | awk '{A[$1]++}END{for(i in A)print i,A[i]}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_gene_number_table.txt

# Calculate the total gene expression in each tile
awk -v OFS='\t' '{print $1, $7}' $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk -F"\t" '{array[$1]+=$2} END { for (i in array) {print i"\t" array[i]}}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_gene_expression_table.txt

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


####################### QC metrics
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

##### Number of input HiFi read pairs
m1=$(grep -w "M\\:\\:mem_process_seqs" $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup_k$bwa_seed_length.log |
cut -d " " -f3 | xargs | tr ' ' + | bc)

##### Number of HiFi-Slide L2R2 passing pear filtering
a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq | cut -d " " -f 1)
m2=$(($a / 4))

a=$(echo "scale=4 ; $m2 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m3=$a"%"

##### Number of HiFi-Slide L2R2 passing length filtering (performed automatically by fastp)
m4=$(grep -w "Number of input reads" $L2R2_GENOME_DIR/L2R2_genome.Log.final.out | cut -d "|" -f2 | sed 's/\t//g')

a=$(echo "scale=4 ; $m4 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m5=$a"%"

##### Number of HiFi-Slide L2R2 uniquely mapped to genome (and to annotated genes)
m6=$(awk '!seen[$1]++' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed | wc -l)

a=$(echo "scale=4 ; $m6 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m7=$a"%"

##### Number of HiFi-Slide L2R1 spatially resolved
m8=$(awk '!seen[$1]++' $L2R1_L1R1_SAM | wc -l)

a=$(echo "scale=4 ; $m8 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m9=$a"%"

##### Number of HiFi read pairs mapped to the genome and spatially resolved
m10=$(awk '!seen[$1]++' $L2R1_MAPPING_DIR/L2R1_L1R1.hifiwrangling0.sort.o | wc -l)

a=$(echo "scale=4 ; $m10 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m11=$a"%"

##### Average number of spots per tile
m12=$(count=0; total=0; for i in $( awk '{ print $2; }' $L2R1_L2R2_INTEGRATE_DIR/tile_spot_number_table.txt );\
do total=$(echo $total+$i | bc ); \
((count++)); done; echo "scale=0; $total / $count" | bc)

m13=$(echo "scale=4 ; $m12 / 10000" | bc | awk '{printf("%.4f",$1)}')

##### Average number of genes per tile
m14=$(count=0; total=0; for i in $( awk '{ print $2; }' $L2R1_L2R2_INTEGRATE_DIR/tile_gene_number_table.txt );\
do total=$(echo $total+$i | bc ); \
((count++)); done; echo "scale=0; $total / $count" | bc)

m15=$(echo "scale=4 ; $m14 / 10000" | bc | awk '{printf("%.4f",$1)}')


##### Print to file
M1="Total number of read pairs"
M2="Number of read pairs passing PEAR filtering"
M3="Percentage of read pairs passing PEAR filtering"
M4="Number of read pairs passing PEAR and FASTP filtering"
M5="Percentage of read pairs passing PEAR and FASTP filtering"
M6="Number of read pairs genome mapped"
M7="Percentage of read pairs genome mapped"
M8="Number of read pairs spatially resolved"
M9="Percentage of read pairs spatially resolved"
M10="Number of read pairs genome mapped and spatially resolved"
M11="Percentage of read pairs genome mapped and spatially resolved"
M12="Average spots per tile"
M13="Average spots per 10 um^2"
M14="Average genes per tile"
M15="Average genes per 10 um^2"

rm $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
touch $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
# for k in $(1 8 9 3 4 5 6 7 10 11 12 13 14 15); do
for k in $(seq 1 15); do
Mk=M${k}
mk=m${k}
echo -e ${!Mk}'\t'${!mk}>> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
done


echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] QC metrics calculation finished." 
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] QC metrics calculation finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


echo "Data processing started: "$START_DATE
echo "Data processing ended: "$(date)

echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo "Data processing started: "$START_DATE >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo "Data processing ended: "$(date) >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log



