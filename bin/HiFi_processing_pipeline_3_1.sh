################## INPUT PARAMETERS (TO BE UPDATED PER EACH SAMPLE)
OUT_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/IGM
SAMPLE_NAME=HiFi_placenta_1
RUNNING_LABEL=""

BIN_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/bin

N_THREADS=32
BWA_MEMORY=80000 # memory (in Megabytes) to be used for bwa index. It does not seem that useful.

mkdir -p $OUT_DIR/$SAMPLE_NAME

# Flowcell and surface identifiers
flowcell_type="NextSeq" # one of: MiniSeq, NextSeq
flowcell="AAALN5GM5"

if [ "$flowcell_type" == "NextSeq" ]; then
surface=$flowcell:1:1
elif [ "$flowcell_type" == "MiniSeq" ]; then
surface=$flowcell:1:
fi

# Directories of the processed data
L1_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/$flowcell # spatial barcodes
L2_DIR=$OUT_DIR/$SAMPLE_NAME # HiFi library

# Spatial barcode files
L1R1_FASTQ_BWA_INDEX=$L1_DIR/bwa_index_L1R1/L1R1_dedup # bwa index path and basename
L1R1_DEDUP=$L1_DIR/L1R1_dedup.fasta # first output of surfdedup
L1R1_DUP=$L1_DIR/L1R1_dup.txt # second output of surfdedup

# Raw reads of HiFi Slides sequencing
L2R1_FASTQ=/mnt/extraids/SDSC_NFS/rcalandrelli/Lab_sequencing/IGM/fastq/221018_A00953_0641_BHGHWNDRX2/JP_placenta_hifi_05oct22_S2_L002_R1_001.fastq.gz
L2R2_FASTQ=/mnt/extraids/SDSC_NFS/rcalandrelli/Lab_sequencing/IGM/fastq/221018_A00953_0641_BHGHWNDRX2/JP_placenta_hifi_05oct22_S2_L002_R2_001.fastq.gz

# Annotation file hg38
annotation_gtf_file=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gtf

# Select full gene coordinates only
annotation_gtf_file_genes=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gene.gtf
# awk -v OFS='\t' '$3=="gene"' $annotation_gtf_file > $annotation_gtf_file_genes

### Mapping reference indexes
STAR_INDEX=/dataOS/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/STARindex_withSJ
BOWTIE2_INDEX=/mnt/extraids/SDSC_NFS/linpei/genome/HSATR
BOWTIE2_INDEX_TRANSCRIPT=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation


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

# ### 1
# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.json \
# --trim_front1 30 \
# --disable_quality_filtering \
# --thread 16

# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.log.json \
# --trim_front1 30 \
# --disable_quality_filtering \
# --thread 16

# ### 2
RUNNING_LABEL="fastp_filter"
fastp \
-i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
-o $L2_DIR/L2R2_preprocessing/L2R2.$RUNNING_LABEL.fastq \
-h $L2_DIR/L2R2_preprocessing/L2R2.$RUNNING_LABEL.log.html \
-j $L2_DIR/L2R2_preprocessing/L2R2.$RUNNING_LABEL.log.json \
--trim_poly_g \
--trim_poly_x \
--disable_quality_filtering \
--thread 16

### 3
# RUNNING_LABEL="fastp_filter"
# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.log.json \
# --trim_poly_g \
# --trim_poly_x \
# --low_complexity_filter \
# --complexity_threshold 30 \
# --disable_quality_filtering \
# --thread 16

### 4
# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.json \
# --trim_front1 30 \
# --disable_quality_filtering \
# --thread 16

# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60_2.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60_2.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60_2.log.json \
# --trim_front1 30 \
# --trim_poly_g \
# --trim_poly_x \
# --disable_quality_filtering \
# --thread 16

# ### 5
# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.json \
# --trim_front1 30 \
# --trim_tail1 30 \
# --disable_quality_filtering \
# --thread 16

# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.trim_front_tail_60.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.trim_front_tail_60.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.trim_front_tail_60.log.json \
# --trim_front1 30 \
# --trim_tail1 30 \
# --disable_quality_filtering \
# --thread 16

# ### 6
# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.json \
# --trim_tail1 30 \
# --disable_quality_filtering \
# --thread 16

# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.trim_tail_60.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.trim_tail_60.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.trim_tail_60.log.json \
# --trim_tail1 30 \
# --disable_quality_filtering \
# --thread 16

# rm $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq

echo "[$(date '+%m-%d-%y %H:%M:%S')] Filtering complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Align HiFi R2 reads to genome/genes in order to obtain gene annotation for HiFi read pairs.
echo "[$(date '+%m-%d-%y %H:%M:%S')] Align HiFi-Slide reads R2 (L2R2) to the genome using STAR..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome
# L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq

L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter
L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.fastq

L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter_ct30_nosjdb
L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.fastq

L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_pear_filter
L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq

# L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter_ct30
# L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.fastq

# L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_2
# L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.trim_front_60_2.fastq

# L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_3
# L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.trim_front_tail_60.fastq

# L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_4
# L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.trim_tail_60.fastq

mkdir -p $L2R2_GENOME_DIR
STAR \
--genomeDir $STAR_INDEX \
--readFilesIn $L2R2_FILTER_FASTQ \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outSAMattributes All \
--outFileNamePrefix $L2R2_GENOME_DIR/L2R2_genome. \
--sjdbGTFfile $annotation_gtf_file \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--runThreadN $N_THREADS


### Select uniquely mapped reads
samtools view -@ $N_THREADS -b -h -q 255 \
-o $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
$L2R2_GENOME_DIR/L2R2_genome.Aligned.sortedByCoord.out.bam

### Map uniquely mapped reads over genes. Cannot use featureCounts because we need to keep track of what L2R2 read align to each gene.
bedtools intersect \
-a $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
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

# rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp.bed
# rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_id.txt
# rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_name.txt
# rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_biotype.txt

### Not aligned to genes
bedtools intersect \
-a $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
-b $annotation_gtf_file_genes \
-v -bed | cut -f 4 | awk '{print $0, "\tNA", "\tNA", "\tNA"}' > $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed

### Concatenate and sort
cat $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_1.bed $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed

rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_1.bed
rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed

echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment of HiFi-Slide reads R2 (L2R2) to the genome complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


### Align HiFi R2 reads to the transcriptome in order to obtain gene annotation for HiFi read pairs.

### Creating Bowtie 2 indexes (if not input parameter)
if [ "$BOWTIE2_INDEX_TRANSCRIPT" == "" ]; then
echo "[$(date '+%m-%d-%y %H:%M:%S')] Create Bowtie 2 indexes of transcriptomes..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

for my_transcript in tRNA piRNA miRNA circRNA; do
mkdir -p /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/$my_transcript/bowtie2_index

bowtie2-build \
--threads $N_THREADS \
--quiet \
/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/$my_transcript/*fa* \
/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/$my_transcript/bowtie2_index/$my_transcript
done
echo "[$(date '+%m-%d-%y %H:%M:%S')] Bowtie 2 index creation complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
fi


# Dummy fastq for testing purposes
# head -40000 $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq > $L2R2_TRANSCRIPTOME_DIR/temp_L2R2.trim_front_60.fastq

### Mapping
echo "[$(date '+%m-%d-%y %H:%M:%S')] Align HiFi-Slide reads R2 (L2R2) to the transcriptomes using Bowtie 2..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

L2R2_TRANSCRIPTOME_DIR=$L2_DIR/L2R2_mapping/transcriptome_fastp_filter
L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.fastq

# L2R2_TRANSCRIPTOME_DIR=$L2_DIR/L2R2_mapping/transcriptome_fastp_filter_ct30
# L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.fastq

mkdir -p $L2R2_TRANSCRIPTOME_DIR

for my_transcript in tRNA piRNA miRNA circRNA; do
mkdir -p $L2R2_TRANSCRIPTOME_DIR/$my_transcript

bowtie2 \
-x $BOWTIE2_INDEX_TRANSCRIPT/$my_transcript/bowtie2_index/$my_transcript \
-U $L2R2_FILTER_FASTQ \
-S $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_mapped.sam" \
--no-unal --threads $N_THREADS --local 2> $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_mapped.log"

# Select uniquely mapped reads. Inverse grep (-v) of reads with auxiliary tag XS, meaning reads that have other valid mappings. This gives exactly the number of reads "aligned exactly 1 time" in the bowtie2 log file.
samtools view $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_mapped.sam" | grep -v "XS:i:" | cut -f 1,3 | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sort.txt"

# rm $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_mapped.sam"

done

echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment of HiFi-Slide reads R2 (L2R2) to the transcriptomes complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide library 2 finished."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide library 2 finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log




#################### LIBRARY 2 R1
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L2R1..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L2R1..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


###### Subsetting barcodes only within tiles under ROI and create specific BWA index
# mkdir -p $L2_DIR/L1R1_ROI
# ROI_tiles=$L2_DIR/L1R1_ROI/ROI_tiles.txt

# less $L1_DIR/L1R1_dedup.fasta | cut -f5 -d ":" | grep -f $ROI_tiles -n | awk -F: '{print $1"\n"$1+1}' > $L2_DIR/L1R1_ROI/temp_lines.txt

# awk 'NR==FNR{data[$1]; next}FNR in data' $L2_DIR/L1R1_ROI/temp_lines.txt $L1_DIR/L1R1_dedup.fasta > $L2_DIR/L1R1_ROI/L1R1_dedup_ROI.fasta

# mkdir -p $L2_DIR/L1R1_ROI/bwa_index_L1R1_ROI
# bwa index \
# -p $L2_DIR/L1R1_ROI/bwa_index_L1R1_ROI/L1R1_dedup_ROI \
# $L2_DIR/L1R1_ROI/L1R1_dedup_ROI.fasta

# L1R1_FASTQ_BWA_INDEX_ROI=$L2_DIR/L1R1_ROI/bwa_index_L1R1_ROI/L1R1_dedup_ROI


###### Align HiFi R1 reads (L2R1) to spatial barcodes (L1R1) in order to obtain spatial coordinates for HiFi read pairs.

echo "[$(date '+%m-%d-%y %H:%M:%S')] Start aligning HiFi-Slide R1 reads (L2R1) to spatial barcodes (L1R1)..." 
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start aligning HiFi-Slide R1 reads (L2R1) to spatial barcodes (L1R1)..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Alignment for all the barcodes
L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping
mkdir -p $L2R1_MAPPING_DIR

# Calculate the minimum base match (-k) as 50% of the length of the spatial barcodes
temp=$(less $L1R1_DEDUP | head -2 | sed -n '2p')
min_base_match=$(echo "scale=4 ; ${#temp} * 0.5" | bc | awk '{printf("%.0f",$1)}')

bwa mem \
-a \
-k $min_base_match \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX $L2R1_FASTQ > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam 2>$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.log


### Alignment for the barcodes under ROI
# L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping_ROI
# mkdir -p $L2R1_MAPPING_DIR

# # Calculate the minimum base match (-k) as 95% of the length of the spatial barcodes
# temp=$(less $L2_DIR/L1R1_ROI/L1R1_dedup_ROI.fasta | head -2 | sed -n '2p')
# min_base_match=$(echo "scale=4 ; ${#temp} * 0.95" | bc | awk '{printf("%.0f",$1)}')

# bwa mem \
# -a \
# -k $min_base_match \
# -t $N_THREADS \
# $L1R1_FASTQ_BWA_INDEX_ROI $L2R1_FASTQ > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam 2>$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.log



### Remove header (not useful and only occupies storage)
grep -v '^@' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.sam
rm $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam

echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment done."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment done." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Select flowcell surface if flowcell type is MiniSeq
if [ "$flowcell_type" == "MiniSeq" ]; then
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selecting the correct flowcell surface..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# Select 0 and 256 flags
awk -F"\t" '$2 == "0" || $2 == "256" { print $0 }' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam

# Select reads coming from surface 1 and surface 2
grep $surface"1" $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface1.sam
grep $surface"2" $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface2.sam

rm $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam

# Choose which surface the tissue is based on the number of mapped reads
n_reads_surface_1=$(wc -l $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface1.sam | cut -d " " -f 1)
n_reads_surface_2=$(wc -l $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface2.sam | cut -d " " -f 1)

if [ $n_reads_surface_1 > $n_reads_surface_2 ]; then
L2R1_L1R1_SAM=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface1.sam
rm $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface2.sam
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selection of the surface done. Selected surface 1." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
else
L2R1_L1R1_SAM=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface2.sam
$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.surface1.sam
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selection of the surface done. Selected surface 2." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
fi

elif [ "$flowcell_type" == "NextSeq" ]; then

L2R1_L1R1_SAM=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.sam

fi


### Filter SAM file to select only HiFi-Slide reads mapped to genome/transcriptome (samf: "SAM filter" custom format)
awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed $L2R1_L1R1_SAM > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.genome.samf

for my_transcript in tRNA piRNA miRNA circRNA; do
size=$(stat -c %s $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sort.txt")

if [ $size != 0 ]; then
awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sort.txt" $L2R1_L1R1_SAM > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.$my_transcript.samf
fi
done

# Concatenate the samf files to obtain a final SAM file
L2R1_L1R1_SAM_FILTER=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.filter.sam
cat $L2R1_MAPPING_DIR/*samf > $L2R1_L1R1_SAM_FILTER

# wc -l $L2R1_L1R1_SAM_FILTER
# cut -f1 $L2R1_L1R1_SAM_FILTER | sort --parallel=$N_THREADS | uniq | wc -l

# rm $L2R1_MAPPING_DIR/*samf
# rm $L2R1_L1R1_SAM

### Select HiFi-Slide R1 reads with highest alignment score
echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score (hifislida.pl)..."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score (hifislida.pl)..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

$BIN_DIR/hifislida.pl $L2R1_L1R1_SAM_FILTER > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida.o 2>$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida.e

echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score (hifislida.pl) complete."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Select HiFi-Slide R1 reads with highest alignment score (hifislida.pl) complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Match HiFi-Slide read pairs with spatial location (header already included)
echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location (hifislida3.pl)..."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location (hifislida3.pl)..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

$BIN_DIR/hifislida3.pl \
$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida.o \
$L1_DIR/L1R1_dup.txt | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort.o

# $BIN_DIR/hifislida3_noDup.pl \
# $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida.o | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort_noDup.o

echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location (hifislida3.pl) complete."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location (hifislida3.pl) complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide L2R1 complete."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide L2R1 complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log




########## Integrate spatial coordinates and gene expression information
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate
mkdir -p $L2R1_L2R2_INTEGRATE_DIR


k_value=k80
L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping/pipeline_3/$k_value/pipeline_3_fastp_filter_ct30
L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter_ct30
L2R2_TRANSCRIPTOME_DIR=$L2_DIR/L2R2_mapping/transcriptome_fastp_filter_ct30
L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate/pipeline_3_fastp_filter_ct30/$k_value/noDup
mkdir -p $L2R1_L2R2_INTEGRATE_DIR

k_value=k95
L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping/pipeline_3/$k_value
L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter_ct30
L2R2_TRANSCRIPTOME_DIR=$L2_DIR/L2R2_mapping/transcriptome_fastp_filter_ct30
L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate/pipeline_3_fastp_filter_ct30/$k_value
mkdir -p $L2R1_L2R2_INTEGRATE_DIR

k_value=k95
L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping_ROI
L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter_ct30
L2R2_TRANSCRIPTOME_DIR=$L2_DIR/L2R2_mapping/transcriptome_fastp_filter_ct30
L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate_ROI
mkdir -p $L2R1_L2R2_INTEGRATE_DIR

### Genome
HiFi_L2R2_genome=$L2R2_GENOME_DIR/HiFi_L2R2_genome.sort.bed

join -1 1 -2 1 -t $'\t' $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort_noDup.o $HiFi_L2R2_genome | cut -f 2,3,4,5,6,7,8 > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt

# rm $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt

# Count how many spots in each tile
awk '{A[$1]++}END{for(i in A)print i,A[i]}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt > $L2R1_L2R2_INTEGRATE_DIR/tile_spot_number_table.txt

# Calculate total gene expression level in each tile
awk '{a[$1]+=$7}END{for(i in a) print i,a[i]}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt > $L2R1_L2R2_INTEGRATE_DIR/tile_gene_expression_table.txt


### TEST
# awk -F"\t" -v OFS='\t' '{ print $1, $2, $3, 1/$4, $5, $6, $7 }' HiFi_L2R2_genome_spatial.txt > HiFi_L2R2_genome_spatial_temp.txt

# awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=$4} END { for (i in array) {print i"\t" array[i]}}' HiFi_L2R2_genome_spatial_temp.txt > HiFi_L2R2_genome_spatial.final2.txt

# cut -f 2,3,4,5,9,10,11 HiFi_L2R2_genome_spatial.txt | tail -n +2 > HiFi_L2R2_genome_spatial2.txt

# tail -n +2 HiFi_L2R2_genome_spatial2.txt > HiFi_L2R2_genome_spatial.txt
# awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' HiFi_L2R2_genome_spatial.txt > HiFi_L2R2_genome_spatial.final.txt


# less HiFi_L2R2_genome_spatial.txt | head -10 > temp.txt
# cut -f 2,3,4,5,9,10,11 temp.txt > temp1.txt
# awk -F"\t" -v OFS='\t' '{ print $1, $2, $3, 1/$4, $5, $6, $7 }' temp1.txt > temp2.txt

# awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=$4} END { for (i in array) {print i"\t" array[i]}}' temp2.txt > temp_agg.txt

# awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' temp1.txt > temp_agg2.txt



### Transcriptome

for my_transcript in tRNA piRNA miRNA circRNA; do

size=$(stat -c %s $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sort.txt")

if [ $size != 0 ]; then # only if there are uniquely mapped reads otherwise do nothing

### COMMENTED FOR UPDATING FILES (THIS IS THE CODE TO BE KEPT)
mkdir -p $L2R1_L2R2_INTEGRATE_DIR/$my_transcript

join -1 1 -2 1 -t $'\t' $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort_noDup.o $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sort.txt" | cut -f 2,3,4,5,6 > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt"

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt"

# Count how many spots in each tile
awk '{A[$1]++}END{for(i in A)print i,A[i]}' $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt" > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/tile_spot_number_table.txt

# Calculate total gene expression level in each tile
awk '{a[$1]+=$7}END{for(i in a) print i,a[i]}' $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt" > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/tile_gene_expression_table.txt

# rm $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt"

### Updating files created with previous pipeline
# cut -f 2,3,4,5,6 $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" | tail -n +2 > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial2.txt"

# awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial2.txt" > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt"

fi

done

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log




###################### Select spots in ROI after visual inspection of the image
# less $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt | head -10 > $L2R1_L2R2_INTEGRATE_DIR/temp.txt

# ROI_TILES=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/result/placenta/HiFi_placenta_1/ROI_tiles.txt



# awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.ROI.txt

# for my_transcript in tRNA piRNA miRNA circRNA; do

# size=$(stat -c %s $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt")

# if [ $size != 0 ]; then # only if there are uniquely mapped reads otherwise do nothing

# awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt" > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.ROI.txt"

# fi
# done






####################### QC metrics (TO BE UPDATED)
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

##### Number of input HiFi read pairs
m1=$(grep -w "M\\:\\:mem_process_seqs" $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.log |
cut -d " " -f3 | xargs | tr ' ' + | bc)

##### Number of HiFi-Slide L2R2 passing pear filtering
a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq | cut -d " " -f 1)
m2=$(($a / 4))

a=$(echo "scale=4 ; $m2 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m3=$a"%"

##### Number of HiFi-Slide L2R2 passing length and low complexity filtering (performed automatically by fastp)
# a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq | cut -d " " -f 1)
# echo $(($a / 4))
# To make it faster, we can just count the number of input reads in the log of the STAR aligner
m4=$(grep -w "Number of input reads" $L2R2_GENOME_DIR/L2R2_genome.Log.final.out |cut -d "|" -f2 | sed 's/\t//g')

a=$(echo "scale=4 ; $m4 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m5=$a"%"

##### Number of HiFi-Slide L2R2 uniquely mapped to genome and to annotated genes
m6=$(cut -f4 $L2R2_GENOME_DIR/HiFi_L2R2_genome.sort.bed | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m6 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m7=$a"%"

##### Print to file
M1="Number of HiFi-Slide read pairs"
M2="Number of HiFi-Slide read pairs passing PEAR filtering"
M3="Percentage of HiFi-Slide read pairs passing PEAR filtering"
M4="Number of HiFi-Slide read pairs passing PEAR and FASTP filtering"
M5="Percentage of HiFi-Slide read pairs passing PEAR and FASTP (length) filtering"
M6="Number of HiFi-Slide read pairs genome mapped (uniquely aligned to annotated genes)"
M7="Percentage of HiFi-Slide read pairs genome mapped"

rm $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
touch $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
for k in $(seq 1 7); do
Mk=M${k}
mk=m${k}
echo -e ${!Mk}'\t'${!mk}>> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
done

##### Number of HiFi-Slide L2R2 uniquely mapped to transcriptome
m9=$m6 # to store number of read pairs uniquely mapped to genome and transcriptome

for my_transcript in tRNA piRNA miRNA circRNA; do

size=$(stat -c %s $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt")

if [ $size != 0 ]; then

##### Number of uniquely mapped reads to the transcriptome
m8a=$(grep -w "aligned exactly 1 time" $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_mapped.log" | cut -d " " -f5)

a=$(echo "scale=4 ; $m8a / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m8b=$a"%"

m9=$(($m9 + $m8a)) # update value

elif [ $size == 0 ]; then
m8a=0
m8b="0%"

fi

M8a=$my_transcript" - Number of HiFi-Slide read pairs transcriptome mapped (uniquely aligned to transcripts)"
M8b=$my_transcript" - Percentage of HiFi-Slide read pairs transcriptome mapped"

for k in a b; do
Mk=M8${k}
mk=m8${k}
echo -e ${!Mk}'\t'${!mk}>> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
done

done

a=$(echo "scale=4 ; $m9 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m10=$a"%"

##### Number of HiFi read pairs aligned to genes/transcripts and spatially resolved
m11=$(wc -l $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort_noDup.o | cut -d " " -f 1)

a=$(echo "scale=4 ; $m11 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m12=$a"%"

M9="Number of HiFi-Slide read pairs genome and transcriptome mapped"
M10="Percentage of HiFi-Slide read pairs genome and transcriptome mapped"
M11="Number of HiFi-Slide read pairs genome and transcriptome mapped and spatially resolved"
M12="Percentage of HiFi-Slide read pairs genome and transcriptome mapped and spatially resolved"

for k in $(seq 9 12); do
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




















