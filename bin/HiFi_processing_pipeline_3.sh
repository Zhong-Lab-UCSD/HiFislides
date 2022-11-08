################## INPUT PARAMETERS (TO BE UPDATED PER EACH SAMPLE)
OUT_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/MiniSeq
SAMPLE_NAME=data_26_05aug22plantrun2

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
L2R1_FASTQ=/mnt/extraids/SDSC_NFS/linpei/hifi/data_26_05aug22plantrun2/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
L2R2_FASTQ=/mnt/extraids/SDSC_NFS/linpei/hifi/data_26_05aug22plantrun2/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz

# Annotation file hg38
annotation_gtf_file=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gtf

### Mapping reference indexes
STAR_INDEX=/dataOS/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/STARindex_withSJ
BOWTIE2_INDEX=/mnt/extraids/SDSC_NFS/linpei/genome/HSATR
BOWTIE2_INDEX_TRANSCRIPT=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation


################## PROCESSING
touch $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

START_DATE=$(date) # start processing date
echo "Processing of "$SAMPLE_NAME
echo "Processing of "$SAMPLE_NAME >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# Select full gene coordinates only
awk -v OFS='\t' '$3=="gene"' $annotation_gtf_file > /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gene.gtf


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
# fastp \
# -i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
# -o $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.fastq \
# -h $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.log.html \
# -j $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.log.json \
# --trim_poly_g \
# --trim_poly_x \
# --disable_quality_filtering \
# --thread 16

### 3
fastp \
-i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
-o $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.fastq \
-h $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.log.html \
-j $L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.log.json \
--trim_poly_g \
--trim_poly_x \
--low_complexity_filter \
--complexity_threshold 30 \
--disable_quality_filtering \
--thread 16

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

# L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter
# L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.fastq

L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter_ct30
L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.fastq

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


# L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter_ct30_066
# L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.fastq
# mkdir -p $L2R2_GENOME_DIR
# STAR \
# --genomeDir $STAR_INDEX \
# --readFilesIn $L2R2_FILTER_FASTQ \
# --outSAMtype BAM SortedByCoordinate \
# --outReadsUnmapped Fastx \
# --outSAMattributes All \
# --outFileNamePrefix $L2R2_GENOME_DIR/L2R2_genome. \
# --sjdbGTFfile $annotation_gtf_file \
# --runThreadN $N_THREADS



### Select uniquely mapped reads
samtools view -@ $N_THREADS -b -h -q 255 \
-o $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
$L2R2_GENOME_DIR/L2R2_genome.Aligned.sortedByCoord.out.bam

### Map uniquely mapped reads over genes. Cannot use featureCounts because we need to keep track of what L2R2 read align to each gene.
bedtools intersect \
-a $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
-b $annotation_gtf_file \
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
$L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_biotype.txt | cut -f 1,3,4,5 | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R2_GENOME_DIR/HiFi_L2R2_genome.sort.bed

rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp.bed
rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_id.txt
rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_name.txt
rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_gene_biotype.txt

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

# L2R2_TRANSCRIPTOME_DIR=$L2_DIR/L2R2_mapping/transcriptome_fastp_filter
# L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter.fastq

L2R2_TRANSCRIPTOME_DIR=$L2_DIR/L2R2_mapping/transcriptome_fastp_filter_ct30
L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.fastp_filter_ct30.fastq

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

### Align HiFi R1 reads (L2R1) to spatial barcodes (L1R1) in order to obtain spatial coordinates for HiFi read pairs.

# Alignment
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start aligning HiFi-Slide R1 reads (L2R1) to spatial barcodes (L1R1)..." 
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start aligning HiFi-Slide R1 reads (L2R1) to spatial barcodes (L1R1)..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping
mkdir -p $L2R1_MAPPING_DIR

# Calculate the minimum base match (-k) as 80% of the length of the spatial barcodes
temp=$(less $L1R1_DEDUP | head -2 | sed -n '2p')
min_base_match=$(echo "scale=4 ; ${#temp} * 0.8" | bc | awk '{printf("%.0f",$1)}')

bwa mem \
-a \
-k $min_base_match \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX $L2R1_FASTQ > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam 2>$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.log

# Remove header (not useful and only occupies storage)
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
awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome.sort.bed $L2R1_L1R1_SAM > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.genome.samf

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

### Rank tiles by the number of HiFi-Slide read pairs
# echo "[$(date '+%m-%d-%y %H:%M:%S')] Rank tiles by the number of HiFi-Slide read pairs (hifislida2.pl)..."
# echo "[$(date '+%m-%d-%y %H:%M:%S')] Rank tiles by the number of HiFi-Slide read pairs (hifislida2.pl)..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# $BIN_DIR/hifislida2.pl \
# $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida.o \
# $L2R1_L1R1_SAM_FILTER > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida2.o 2>$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida2.e

# echo "[$(date '+%m-%d-%y %H:%M:%S')] Rank tiles by the number of HiFi-Slide read pairs (hifislida2.pl) complete."
# echo "[$(date '+%m-%d-%y %H:%M:%S')] Rank tiles by the number of HiFi-Slide read pairs (hifislida2.pl) complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# ### Select tiles under ROI
# echo "[$(date '+%m-%d-%y %H:%M:%S')] Select tiles under ROI..."
# echo "[$(date '+%m-%d-%y %H:%M:%S')] Select tiles under ROI..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# if [ "$flowcell_type" == "MiniSeq" ]; then
# if [ $n_reads_surface_1 > $n_reads_surface_2 ]; then
# mySurf=1
# else
# mySurf=2
# fi
# elif [ "$flowcell_type" == "NextSeq" ]; then
# mySurf=1
# fi

# $BIN_DIR/select_tiles_in_ROI.r \
# -i $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida2.o \
# -o $L2R1_MAPPING_DIR/ROI_tile_IDs.txt \
# -f $flowcell_type \
# --surface $mySurf \
# --max_size_ROI $max_size_ROI \
# --min_size_ROI $min_size_ROI \
# --p_value 0.1

# echo "[$(date '+%m-%d-%y %H:%M:%S')] ROI selection complete."
# echo "[$(date '+%m-%d-%y %H:%M:%S')] ROI selection complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Match HiFi-Slide read pairs with spatial location (header already included)
echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location (hifislida3.pl)..."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location (hifislida3.pl)..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# cut -f2 $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida2.o | sort --parallel=$N_THREADS > $L2R1_MAPPING_DIR/ALL_tile_IDs.txt 

$BIN_DIR/hifislida3.pl \
$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida.o \
$L1_DIR/L1R1_dup.txt | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort.o


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


### TEMPORARY
L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping/pipeline_3/pipeline_3_fastp_filter_ct30
L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter_ct30
L2R2_TRANSCRIPTOME_DIR=$L2_DIR/L2R2_mapping/transcriptome_fastp_filter_ct30
L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate/pipeline_3_fastp_filter_ct30

L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping/pipeline_3/pipeline_3_fastp_filter
L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter
L2R2_TRANSCRIPTOME_DIR=$L2_DIR/L2R2_mapping/transcriptome_fastp_filter
L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate/pipeline_3_fastp_filter


### Genome
HiFi_L2R2_genome=$L2R2_GENOME_DIR/HiFi_L2R2_genome.sort.bed

join -1 1 -2 1 -t $'\t' $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort.o $HiFi_L2R2_genome | cut -f 2,3,4,5,6,7,8 > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt

# rm $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt

# Count how many spots in each tile
awk '{A[$1]++}END{for(i in A)print i,A[i]}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt > $L2R1_L2R2_INTEGRATE_DIR/tile_spot_number_table.txt


### TEST
awk -F"\t" -v OFS='\t' '{ print $1, $2, $3, 1/$4, $5, $6, $7 }' HiFi_L2R2_genome_spatial.txt > HiFi_L2R2_genome_spatial_temp.txt

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=$4} END { for (i in array) {print i"\t" array[i]}}' HiFi_L2R2_genome_spatial_temp.txt > HiFi_L2R2_genome_spatial.final2.txt

cut -f 2,3,4,5,9,10,11 HiFi_L2R2_genome_spatial.txt | tail -n +2 > HiFi_L2R2_genome_spatial2.txt

tail -n +2 HiFi_L2R2_genome_spatial2.txt > HiFi_L2R2_genome_spatial.txt
awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' HiFi_L2R2_genome_spatial.txt > HiFi_L2R2_genome_spatial.final.txt


less HiFi_L2R2_genome_spatial.txt | head -10 > temp.txt
cut -f 2,3,4,5,9,10,11 temp.txt > temp1.txt
awk -F"\t" -v OFS='\t' '{ print $1, $2, $3, 1/$4, $5, $6, $7 }' temp1.txt > temp2.txt

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=$4} END { for (i in array) {print i"\t" array[i]}}' temp2.txt > temp_agg.txt

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' temp1.txt > temp_agg2.txt



### Transcriptome (FILES TO BE UPDATED)

for my_transcript in tRNA piRNA miRNA circRNA; do

size=$(stat -c %s $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sort.txt")

if [ $size != 0 ]; then # only if there are uniquely mapped reads otherwise do nothing

### COMMENTED FOR UPDATING FILES (THIS IS THE CODE TO BE KEPT)
mkdir -p $L2R1_L2R2_INTEGRATE_DIR/$my_transcript

join -1 1 -2 1 -t $'\t' $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort.o $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sort.txt" | cut -f 2,3,4,5,6 > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt"

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt"

# rm $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt"


### Updating files created with previous pipeline
# cut -f 2,3,4,5,6 $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" | tail -n +2 > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial2.txt"

# awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial2.txt" > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt"

fi

done

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log




###################### Select spots in ROI after visual inspection of the image
less $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt | head -10 > $L2R1_L2R2_INTEGRATE_DIR/temp.txt

ROI_TILES=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/result/placenta/HiFi_placenta_1/ROI_tiles.txt



awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.ROI.txt

for my_transcript in tRNA piRNA miRNA circRNA; do

size=$(stat -c %s $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt")

if [ $size != 0 ]; then # only if there are uniquely mapped reads otherwise do nothing

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.txt" > $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.final.ROI.txt"

fi
done






####################### QC metrics (TO BE UPDATED)
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

##### Number of input HiFi read pairs
m4=$(grep -w "M\\:\\:mem_process_seqs" $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.log |
cut -d " " -f3 | xargs | tr ' ' + | bc)

##### Number of HiFi read pairs aligned to spatial barcodes (number of unique HiFi read IDs in column 1 of L2R1_L1R1_dedup.hifislida.o)
m5=$(awk -F '\t' '{print $1}' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida.o | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m5 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m6=$a"%"


##### Number of HiFi read pairs aligned to unique spatial barcodes (number of HiFi read IDs in column 1 of L2R1_L1R1_dedup.hifislida.o with column 3 and 4 equal to 1)
m7=$(cut -f3 $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida2.o | xargs | tr ' ' + | bc) # same as below but faster
# awk -F "\t" '$3 == 1 && $4 == 1 { print $1 }' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup_1_1.hifislida.o | wc -l

a=$(echo "scale=4 ; $m7 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m8=$a"%"

##### Number of HiFi read pairs aligned to unique spatial barcodes and under ROI
# awk 'BEGIN{FS=OFS="\t"} {gsub(/T/, "", $2)} 1' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup_1_1.hifislida2.o | head -10 # to remove the T (not used anymore)
m9=$(awk -F "\t" 'NR==FNR{a[$1]; next} FNR==0 || $2 in a' $L2R1_MAPPING_DIR/ROI_tile_IDs.txt $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida2.o | cut -f3 | xargs | tr ' ' + | bc)

a=$(echo "scale=4 ; $m9 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m10=$a"%"


##### Number of HiFi read pairs under ROI (all spatially resolved, i.e. not only aligned to unique spatial barcodes)
m11=$(awk -F '\t' 'NR>1 {print $1}' $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.o | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m11 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m12=$a"%"


##### Number of HiFi-Slide L2R2 passing pear filtering
a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2_pear.unassembled.reverse.fastq | cut -d " " -f 1)
m13=$(($a / 4))

a=$(echo "scale=4 ; $m13 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m14=$a"%"

##### Number of HiFi-Slide L2R2 passing length filtering (performed automatically by fastp)
# a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq | cut -d " " -f 1)
# echo $(($a / 4))
# To make it faster, we can just count the number of input reads in the log of the STAR aligner
m15=$(grep -w "Number of input reads" $L2R2_GENOME_DIR/L2R2_genome.Log.final.out |cut -d "|" -f2 | sed 's/\t//g')

a=$(echo "scale=4 ; $m15 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m16=$a"%"


########## GENOME

##### Number of HiFi-Slide L2R2 uniquely mapped to genome
# grep -w "Uniquely mapped reads number" $L2R2_GENOME_DIR/L2R2_genome.Log.final.out | cut -d "|" -f2 | sed 's/\t//g'

##### Number of HiFi-Slide L2R2 uniquely mapped to genome and to annotated genes
m17=$(cut -f4 $L2R2_GENOME_DIR/HiFi_L2R2_genome.bed | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m17 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m18=$a"%"

##### Number of HiFi-Slide read pairs genome mapped and spatially resolved
m19=$(awk 'FNR==NR {a[$1]; next} FNR> 0 && $4 in a' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida.o $L2R2_GENOME_DIR/HiFi_L2R2_genome.bed | cut -f4 | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m19 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m20=$a"%"


##### Number of HiFi reads spatially resolved, under ROI and aligned to genome
m21=$(awk -F '\t' 'NR>1 {print $1}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m21 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m22=$a"%"

# Number of genes 
m23=$(awk -F '\t' 'NR>1 {print $9}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt | sort --parallel=$N_THREADS | uniq | wc -l)

##### Average number of HiFi-Slide read pairs genome mapped and spatially resolved per tile
n_tiles=$(wc -l $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida2.o | cut -d " " -f 1)
m24=$(echo "scale=4 ; $m19 / $n_tiles" | bc | awk '{printf("%.2f",$1)}')

##### Average number of HiFi-Slide read pairs genome mapped and spatially resolved per tile under ROI
n_tiles_under_ROI=$(wc -l $L2R1_MAPPING_DIR/ROI_tile_IDs.txt | cut -d " " -f 1)
m25=$(echo "scale=4 ; $m21 / $n_tiles_under_ROI" | bc | awk '{printf("%.2f",$1)}')


#### Labels
M4="Number of HiFi-Slide read pairs"
M5="Number of HiFi-Slide read pairs spatially resolved (aligned to spatial barcodes)"
M6="Percentage of HiFi-Slide read pairs spatially resolved"
M7="Number of HiFi-Slide read pairs univocally spatially resolved (aligned to unique spatial barcodes)"
M8="Percentage of HiFi-Slide read pairs univocally spatially resolved"
M9="Number of HiFi-Slide read pairs univocally spatially resolved and under ROI"
M10="Percentage of HiFi-Slide read pairs univocally spatially resolved and under ROI"
M11="Number of HiFi-Slide read pairs spatially resolved (aligned to spatial barcodes) and under ROI"
M12="Percentage of HiFi-Slide read pairs spatially resolved (aligned to spatial barcodes) and under ROI"
M13="Number of HiFi-Slide read pairs passing PEAR filtering"
M14="Percentage of HiFi-Slide read pairs passing PEAR filtering"
M15="Number of HiFi-Slide read pairs passing PEAR and FASTP (length) filtering"
M16="Percentage of HiFi-Slide read pairs passing PEAR and FASTP (length) filtering"
M17="Number of HiFi-Slide read pairs genome mapped (uniquely aligned to annotated genes)"
M18="Percentage of HiFi-Slide read pairs genome mapped"
M19="Number of HiFi-Slide read pairs genome mapped and spatially resolved"
M20="Percentage of HiFi-Slide read pairs genome mapped and spatially resolved"
M21="Number of HiFi-Slide read pairs genome mapped and under ROI"
M22="Percentage of HiFi-Slide read pairs genome mapped and under ROI"
M23="Number of annotated genes aligned to HiFi-Slide read pairs"
M24="Average number of HiFi-Slide read pairs genome mapped and spatially resolved per tile"
M25="Average number of HiFi-Slide read pairs genome mapped per tile under ROI"


rm $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
touch $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
for k in $(seq 4 25); do
Mk=M${k}
mk=m${k}
echo -e ${!Mk}'\t'${!mk}>> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
done


########## TRANSCRIPTOME
m27=$m19 # to store number of read pairs genome and transcriptome mapped and spatially resolved
m29=$m21 # to store number of read pairs genome and transcriptome mapped and undr ROI

for my_transcript in tRNA piRNA miRNA circRNA; do

size=$(stat -c %s $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt")

if [ $size != 0 ]; then

##### Number of uniquely mapped reads to the transcriptome
m26a=$(grep -w "aligned exactly 1 time" $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_mapped.log" | cut -d " " -f5)

a=$(echo "scale=4 ; $m26a / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m26b=$a"%"

##### Number of HiFi-Slide read pairs transcriptome mapped and spatially resolved
m26c=$(awk 'FNR==NR {a[$1]; next} FNR> 0 && $1 in a' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida.o $L2R2_TRANSCRIPTOME_DIR/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt" | cut -f1 | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m26c / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m26d=$a"%"

m27=$(($m27 + $m26c)) # update value

##### Number of HiFi reads spatially resolved, under ROI and aligned to transcriptome
m26e=$(awk -F '\t' 'NR>1 {print $1}' $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m26c / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m26f=$a"%"

m29=$(($m29 + $m26e)) # update value

# Number of transcripts 
m26g=$(awk -F '\t' 'NR>1 {print $6}' $L2R1_L2R2_INTEGRATE_DIR/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" | sort --parallel=$N_THREADS | uniq | wc -l)

##### Average number of HiFi-Slide read pairs transcriptome mapped and spatially resolved per tile
n_tiles=$(wc -l $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.hifislida2.o | cut -d " " -f 1)
m26h=$(echo "scale=4 ; $m26c / $n_tiles" | bc | awk '{printf("%.2f",$1)}')

##### Average number of HiFi-Slide read pairs transcriptome mapped and spatially resolved per tile under ROI
n_tiles_under_ROI=$(wc -l $L2R1_MAPPING_DIR/ROI_tile_IDs.txt | cut -d " " -f 1)
m26i=$(echo "scale=4 ; $m26e / $n_tiles_under_ROI" | bc | awk '{printf("%.2f",$1)}')

elif [ $size == 0 ]; then
m26a=0
m26b="0%"
m26c=0
m26d="0%"
m26e=0
m26f="0%"
m26g=0
m26h=0
m26i=0

fi

M26a=$my_transcript" - Number of HiFi-Slide read pairs transcriptome mapped (uniquely aligned to transcripts)"
M26b=$my_transcript" - Percentage of HiFi-Slide read pairs transcriptome mapped"
M26c=$my_transcript" - Number of HiFi-Slide read pairs transcriptome mapped and spatially resolved"
M26d=$my_transcript" - Percentage of HiFi-Slide read pairs transcriptome mapped and spatially resolved"
M26e=$my_transcript" - Number of HiFi-Slide read pairs transcriptome mapped and under ROI"
M26f=$my_transcript" - Percentage of HiFi-Slide read pairs transcriptome mapped and under ROI"
M26g=$my_transcript" - Number of transcripts aligned to HiFi-Slide read pairs"
M26h=$my_transcript" - Average number of HiFi-Slide read pairs transcriptome mapped and spatially resolved per tile"
M26i=$my_transcript" - Average number of HiFi-Slide read pairs transcriptome mapped per tile under ROI"

for k in a b c d e f g h i; do
Mk=M26${k}
mk=m26${k}
echo -e ${!Mk}'\t'${!mk}>> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
done

done

######## Genome and transcriptome merged metrics
a=$(echo "scale=4 ; $m27 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m28=$a"%"

a=$(echo "scale=4 ; $m29 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m30=$a"%"

m31=$(echo "scale=4 ; $m27 / $n_tiles" | bc | awk '{printf("%.2f",$1)}')
m32=$(echo "scale=4 ; $m29 / $n_tiles_under_ROI" | bc | awk '{printf("%.2f",$1)}')

M27="Number of HiFi-Slide read pairs genome and transcriptome mapped and spatially resolved"
M28="Percentage of HiFi-Slide read pairs genome and transcriptome mapped and spatially resolved"
M29="Number of HiFi-Slide read pairs genome and transcriptome mapped and under ROI"
M30="Percentage of HiFi-Slide read pairs genome and transcriptome mapped and under ROI"
M31="Average number of HiFi-Slide read pairs genome and transcriptome mapped and spatially resolved per tile"
M32="Average number of HiFi-Slide read pairs genome and transcriptome mapped per tile under ROI"

for k in $(seq 27 32); do
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




















