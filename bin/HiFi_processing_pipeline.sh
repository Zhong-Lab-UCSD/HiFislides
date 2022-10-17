################## INPUT PARAMETERS
OUT_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data
SAMPLE_NAME=data14_test

BIN_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/bin

N_THREADS=32
BWA_MEMORY=80000 # memory (in Megabytes) to be used for bwa index. It does not seem that useful.

mkdir -p $OUT_DIR/$SAMPLE_NAME

# Directories of the processed data
L1_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/data14_test/lib1 # spatial barcodes
L2_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/data14_test/lib2 # HiFi library
mkdir -p $L1_DIR
mkdir -p $L2_DIR

# Directories of the raw fastq files for each library. The full path is used here.
L1_FASTQ_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/test_sample/lib1/fastq
L1_FASTQ_BASENAME=MT*_L001_R1_001.fastq.gz
L2_FASTQ_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/test_sample/lib2/fastq

# Raw reads of HiFi Slides sequencing
L2R1_FASTQ=$L2_FASTQ_DIR/Undetermined_S0_L001_R1_001.fastq.gz
L2R2_FASTQ=$L2_FASTQ_DIR/Undetermined_S0_L001_R2_001.fastq.gz

# Flowcell and surface identifiers
flowcell_type="NextSeq" # one of: MiniSeq, NextSeq
flowcell=AAAL33WM5

if [ $flowcell_type == "NextSeq" ]; then
surface=$flowcell:1:1
elif [ $flowcell_type == "MiniSeq" ]; then
surface=$flowcell:1:
fi

# Annotation file hg38. This file can be downloaded without the need of computing it from the GTF file or bedtools intersect can take GTF as input, only genes can be selected from the GTF file.
# annotation_gtf_file=/dataOS/sysbio/Genomes/Homo_sapiens/Ensembl/GRCH38_hg38/Annotation/Genes/Homo_sapiens.GRCh38.84.chr.gtf
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
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# Select full genes only
# awk -v OFS='\t' '$3=="gene"' $annotation_gtf_file > /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/Homo_sapiens.GRCh38.84.chr.gene.gtf
awk -v OFS='\t' '$3=="gene"' $annotation_gtf_file > /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gene.gtf


#################### LIBRARY 1 (spatial barcodes)
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide library 1..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide library 1..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Deduplication of raw reads from the recycled flow cell to extract unique raw reads as spatial barcodes
# g++ surfdedup.cpp -o surfdedup -lz
if [ -f "$L1_DIR/L1R1_dedup.fasta" ] && [ -f "$L1_DIR/L1R1_dup.txt" ]; then
echo "[$(date '+%m-%d-%y %H:%M:%S')] L1R1 reads already processed."
echo "[$(date '+%m-%d-%y %H:%M:%S')] L1R1 reads already processed." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
else
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start deduplication of L1R1 reads..."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start deduplication of L1R1 reads..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
$BIN_DIR/surfdedup $surface $L1_FASTQ_DIR/$L1_FASTQ_BASENAME > $L1_DIR/L1R1_dedup.fasta 2>$L1_DIR/L1R1_dup.txt
echo "[$(date '+%m-%d-%y %H:%M:%S')] Deduplication of L1R1 reads complete."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Deduplication of L1R1 reads complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
fi

# dummy example to make running faster
# $BIN_DIR/surfdedup $surface $L1_FASTQ_DIR/MT080_S1_L001_R1_001.fastq.gz > $L1_DIR/L1R1_dedup.fasta 2>$L1_DIR/L1R1_dup.txt


### Align HiFi R1 reads (L2R1) to spatial barcodes (L1R1) in order to obtain spatial coordinates for HiFi read pairs.

# Create index files for L1R1
if ls $L1_DIR/bwa_index_L1R1/${L1R1_dedup}* > /dev/null 2>&1; then
echo "[$(date '+%m-%d-%y %H:%M:%S')] BWA index for spatial barcodes (L1R1) already existing."
echo "[$(date '+%m-%d-%y %H:%M:%S')] BWA index for spatial barcodes (L1R1) already existing." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
else
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start creating BWA index for spatial barcodes (L1R1)..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
BWA_BLOCK_SIZE=$(($BWA_MEMORY * 1000000 / 8)) # currently not used
mkdir -p $L1_DIR/bwa_index_L1R1
bwa index \
-p $L1_DIR/bwa_index_L1R1/L1R1_dedup \
$L1_DIR/L1R1_dedup.fasta
echo "[$(date '+%m-%d-%y %H:%M:%S')] BWA index creation complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
fi

# Alignment
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start aligning HiFi-Slide R1 reads (L2R1) to spatial barcodes (L1R1)..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
mkdir -p $L2_DIR/L2R1_mapping
bwa mem -a -k 40 -t $N_THREADS $L1_DIR/bwa_index_L1R1/L1R1_dedup $L2R1_FASTQ > $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.sam 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.log
echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment done." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

if [ $flowcell_type == "MiniSeq" ]; then
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selecting the correct flowcell surface..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# Remove header and select 0 and 256 flags
grep -v '^@' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.sam | awk -F"\t" '$2 == "0" || $2 == "256" { print $0 }' > $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.temp.sam

# Select reads coming from surface 1 and surface 2
grep $surface"1" $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.temp.sam > $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.surface1.sam
grep $surface"2" $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.temp.sam > $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.surface2.sam

# Choose which surface the tissue is based on the number of mapped reads
n_reads_surface_1=$(wc -l $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.surface1.sam | cut -d " " -f 1)
n_reads_surface_2=$(wc -l $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.surface2.sam | cut -d " " -f 1)

if [ $n_reads_surface_1 > $n_reads_surface_2 ]; then
L2R1_L1R1_SAM=$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.surface1.sam
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selection of the surface done. Selected surface 1." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
else
L2R1_L1R1_SAM=$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.surface2.sam
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selection of the surface done. Selected surface 2." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
fi

rm $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.temp.sam

elif [ $flowcell_type == "NextSeq" ]; then
L2R1_L1R1_SAM=$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.sam
fi


### Select HiFi-Slide R1 reads with highest alignment score
echo "[$(date '+%m-%d-%y %H:%M:%S')] Parse aligned HiFi-Slide R1 reads (L1R1) and select ROI..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
$BIN_DIR/hifislida.pl $L2R1_L1R1_SAM > $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.o 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.e

### Rank the tiles by number of HiFi-Slide read pairs
$BIN_DIR/hifislida2.pl \
$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.o \
$L2R1_L1R1_SAM > $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o 2>$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.e

### Select tiles under ROI
if [ $flowcell_type == "MiniSeq" ]; then
if [ $n_reads_surface_1 > $n_reads_surface_2 ]; then
mySurf=1
else
mySurf=2
fi
elif [ $flowcell_type == "NextSeq" ]; then
mySurf=1
fi

$BIN_DIR/select_tiles_in_ROI.r \
-i $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o \
-o $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt \
-f $flowcell_type \
--surface $mySurf \
--max_size_ROI 4 \
--min_size_ROI 2 \
--p_value 0.05

echo "[$(date '+%m-%d-%y %H:%M:%S')] Parsing and ROI selection complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Match HiFi-Slide read pairs with spatial location
echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
$BIN_DIR/hifislida3.pl \
$L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.o \
$L2_DIR/L2R1_mapping/ROI_tile_IDs.txt \
$L1_DIR/L1R1_dup.txt > $L2_DIR/L2R1_mapping/temp.hifislida3.o

# Add header
echo -e "HiFi_read_id\ttile_id\tcol\trow\tN" | cat - $L2_DIR/L2R1_mapping/temp.hifislida3.o > $L2_DIR/L2R1_mapping/L2R1_L1R1.hifislida3.o

rm $L2_DIR/L2R1_mapping/temp.hifislida3.o
echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide library 1 complete."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide library 1 complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


#################### LIBRARY 2 (HiFi-Slide read pairs)
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide library 2..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide library 2..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Preprocessing of HiFi R2 reads
mkdir -p $L2_DIR/L2R2_preprocessing

# Find L2R2 reads overlapping L2R1 and filter them out using the software pear
echo "[$(date '+%m-%d-%y %H:%M:%S')] Find L2R2 reads overlapping L2R1 and filter them out using the software PEAR..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
minoverlap=10

pear \
-f $L2R1_FASTQ \
-r $L2R2_FASTQ \
-v $minoverlap \
-j $N_THREADS \
-o $L2_DIR/L2R2_preprocessing/L2R2_pear

# grep @ $L2_DIR/L2R2_preprocessing/L2R2_pear.unassembled.reverse.fastq > $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.names
# zgrep -A 3 -f $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.names $L2R2_FASTQ > $L2_DIR/L2R2_preprocessing/L2R2.pear_filter_names.fastq

# Necessary? To be confirmed!
seqtk seq -r $L2_DIR/L2R2_preprocessing/L2R2_pear.unassembled.reverse.fastq > $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq
echo "[$(date '+%m-%d-%y %H:%M:%S')] PEAR processing complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# Trimming the front 60 bp to remove Illumina adapters (max 16 threads allowed, max 30 bp at the time)
echo "[$(date '+%m-%d-%y %H:%M:%S')] Trimming the front 60 bp of L2R2 reads to remove Illumina adapters..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
fastp \
-i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
-o $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
-h $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.html \
-j $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.log.json \
--trim_front1 30 \
--disable_quality_filtering \
--thread 16

fastp \
-i $L2_DIR/L2R2_preprocessing/L2R2.trim_front_temp.fastq \
-o $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq \
-h $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.log.html \
-j $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.log.json \
--trim_front1 30 \
--disable_quality_filtering \
--thread 16
echo "[$(date '+%m-%d-%y %H:%M:%S')] Trimming complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


### Align HiFi R2 reads to genome/genes in order to obtain gene annotation for HiFi read pairs.
echo "[$(date '+%m-%d-%y %H:%M:%S')] Align HiFi-Slide reads R2 (L2R2) to the genome using STAR..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
mkdir -p $L2_DIR/L2R2_mapping/genome

STAR \
--genomeDir $STAR_INDEX \
--readFilesIn $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outSAMattributes All \
--outFileNamePrefix $L2_DIR/L2R2_mapping/genome/L2R2_genome. \
--sjdbGTFfile $annotation_gtf_file \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--runThreadN $N_THREADS

### Select uniquely mapped reads
samtools view -@ $N_THREADS -b -h -q 255 \
-o $L2_DIR/L2R2_mapping/genome/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
$L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam

# samtools view -@ 32 $L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam | wc -l
# samtools view -@ 32 -q 255 $L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam | wc -l
# samtools view -@ 32 -q 30 $L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam | wc -l


### Map uniquely mapped reads over genes. Cannot use featureCounts because we need to keep track of what L2R2 read align to each gene.
bedtools intersect \
-a $L2_DIR/L2R2_mapping/genome/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
-b /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gene.gtf \
-wb -bed | cut -f 1,2,3,4,21 > $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_temp.bed

cut -f 5 $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_temp.bed |  # extract column 5
cut -f1 -d';' |        # extract column 1 using semi-colon as delimiter
cut -f2 -d' ' |        # extract column 2 using space as delimiter
tr -d '"' > $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_gene_id.txt  # remove double quote

cut -f 5 $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_temp.bed |  # extract column 5
cut -f3 -d';' |        # extract column 3 using semi-colon as delimiter
cut -f3 -d' ' |        # extract column 3 using space as delimiter
tr -d '"' > $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_gene_name.txt  # remove double quote

cut -f 5 $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_temp.bed |  # extract column 5
cut -f2 -d';' |        # extract column 2 (column 5 for Ensembl) using semi-colon as delimiter
cut -f3 -d' ' |        # extract column 3 using space as delimiter
tr -d '"' > $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_gene_biotype.txt  # remove double quote

paste \
$L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_temp.bed \
$L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_gene_id.txt \
$L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_gene_name.txt \
$L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_gene_biotype.txt |
cut -f 1,2,3,4,6,7,8 > $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed
echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment of HiFi-Slide reads R2 (L2R2) to the genome complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# rm $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_temp.bed
# rm $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_gene_id.txt
# rm $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_gene_name.txt
# rm $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_gene_biotype.txt


### Align HiFi R2 reads to the transcriptome in order to obtain gene annotation for HiFi read pairs.
mkdir -p $L2_DIR/L2R2_mapping/transcriptome

### Creating Bowtie 2 indexes (if not input parameter)

if [ $BOWTIE2_INDEX_TRANSCRIPT == "" ]; then
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
# head -40000 $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq > $L2_DIR/L2R2_mapping/transcriptome/temp_L2R2.trim_front_60.fastq

### Mapping
echo "[$(date '+%m-%d-%y %H:%M:%S')] Align HiFi-Slide reads R2 (L2R2) to the transcriptomes using Bowtie 2..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
for my_transcript in tRNA piRNA miRNA circRNA; do
mkdir -p $L2_DIR/L2R2_mapping/transcriptome/$my_transcript

bowtie2 \
-x $BOWTIE2_INDEX_TRANSCRIPT/$my_transcript/bowtie2_index/$my_transcript \
-U $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq \
-S $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_mapped.sam" \
--un $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_unmapped.txt" \
--no-unal --threads $N_THREADS --local 2> $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_mapped.log"

# Select uniquely mapped reads
# Option 1: Inverse grep (-v) of reads with auxiliary tag XS, meaning reads that have other valid mappings. This gives exactly the number of reads "aligned exactly 1 time" in the bowtie2 log file.
samtools view $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_mapped.sam" | grep -v "XS:i:" > $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sam"

# Option 2: using MAPQ value
# samtools view -q 10 \
# -o $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_mapped.mapq10.sam" \
# $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_mapped.sam"

# Extract fields of interest
cut -f 1,3 $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sam" > $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt"

done

echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment of HiFi-Slide reads R2 (L2R2) to the transcriptomes complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide library 2 finished."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide library 2 finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log



########## Integrate spatial coordinates and gene expression information
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

mkdir -p $L2_DIR/L2R1_L2R2_integrate

HiFi_L2R1_spatial=$L2_DIR/L2R1_mapping/L2R1_L1R1.hifislida3.sort.o
#sort -k 1 $L2_DIR/L2R1_mapping/hifislida3.o > $HiFi_L2R1_spatial
cat $L2_DIR/L2R1_mapping/L2R1_L1R1.hifislida3.o | sort -k 1 --parallel=$N_THREADS -S 20G > $HiFi_L2R1_spatial

### Genome
HiFi_L2R2_genome=$L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.sort.bed
#sort -k 4 $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed > $HiFi_L2R2_genome
cat $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed | sort -k 4 --parallel=$N_THREADS -S 20G > $HiFi_L2R2_genome

join -1 1 -2 4 -t $'\t' $HiFi_L2R1_spatial $HiFi_L2R2_genome > $L2_DIR/L2R1_L2R2_integrate/temp_HiFi_L2R2_genome_spatial.txt

# Add header
echo -e "HiFi_read_id\ttile_id\tcol\trow\tN\tHiFi_read_chr\tHiFi_read_start\tHiFi_read_end\tgene_id\tgene_name\tgene_type" | cat - $L2_DIR/L2R1_L2R2_integrate/temp_HiFi_L2R2_genome_spatial.txt > $L2_DIR/L2R1_L2R2_integrate/HiFi_L2R2_genome_spatial.txt

rm $L2_DIR/L2R1_L2R2_integrate/temp_HiFi_L2R2_genome_spatial.txt

### Transcriptome

for my_transcript in tRNA piRNA miRNA circRNA; do

size=$(stat -c %s $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt")

if [ $size != 0 ]; then # only if there are uniquely mapped reads otherwise do nothing

mkdir -p $L2_DIR/L2R1_L2R2_integrate/$my_transcript

cat $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt" | sort -k 1 --parallel=$N_THREADS -S 20G > $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sort.txt"

join -1 1 -2 1 -t $'\t' $HiFi_L2R1_spatial $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.sort.txt" > $L2_DIR/L2R1_L2R2_integrate/$my_transcript/temp_HiFi_L2R2_$my_transcript"_spatial.txt"

# Add header
echo -e "HiFi_read_id\ttile_id\tcol\trow\tN\ttranscript_id" | cat - $L2_DIR/L2R1_L2R2_integrate/$my_transcript/temp_HiFi_L2R2_$my_transcript"_spatial.txt" > $L2_DIR/L2R1_L2R2_integrate/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt"

rm $L2_DIR/L2R1_L2R2_integrate/$my_transcript/temp_HiFi_L2R2_$my_transcript"_spatial.txt"

fi

done

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


####################### QC metrics
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

##### Total number of barcodes (L1R1)
rm $L1_DIR/L1R1_stats.txt
touch $L1_DIR/L1R1_stats.txt

for x in $L1_FASTQ_DIR/$L1_FASTQ_BASENAME; do
echo $x
a=$(unpigz -p $N_THREADS -c $x | wc -l)
out=$(($a / 4))
echo -e $x"\t"$out >> $L1_DIR/L1R1_stats.txt
done

m1=$(awk '{ sum += $2 } END { print sum }' $L1_DIR/L1R1_stats.txt)

##### Number of deduplicated barcodes
a=$(wc -l $L1_DIR/L1R1_dedup.fasta | cut -d " " -f 1)
m2=$(($a / 2))

a=$(echo "scale=4 ; $m2 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m3=$a"%"

##### Number of input HiFi read pairs
m4=$(grep -w "M\\:\\:mem_process_seqs" $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.log |
cut -d " " -f3 | xargs | tr ' ' + | bc)

##### Number of HiFi read pairs aligned to spatial barcodes (number of unique HiFi read IDs in column 1 of L2R1_L1R1_dedup.hifislida.o)
m5=$(awk -F '\t' '{print $1}' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.o | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m5 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m6=$a"%"


##### Number of HiFi read pairs aligned to unique spatial barcodes (number of HiFi read IDs in column 1 of L2R1_L1R1_dedup.hifislida.o with column 3 and 4 equal to 1)
m7=$(cut -f3 $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o | xargs | tr ' ' + | bc) # same as below but faster
# awk -F "\t" '$3 == 1 && $4 == 1 { print $1 }' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup_1_1.hifislida.o | wc -l

a=$(echo "scale=4 ; $m7 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m8=$a"%"

##### Number of HiFi read pairs aligned to unique spatial barcodes and under ROI
# awk 'BEGIN{FS=OFS="\t"} {gsub(/T/, "", $2)} 1' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup_1_1.hifislida2.o | head -10 # to remove the T (not used anymore)
m9=$(awk -F "\t" 'NR==FNR{a[$1]; next} FNR==0 || $2 in a' $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o | cut -f3 | xargs | tr ' ' + | bc)

a=$(echo "scale=4 ; $m9 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m10=$a"%"


##### Number of HiFi read pairs under ROI (all spatially resolved, i.e. not only aligned to unique spatial barcodes)
m11=$(awk -F '\t' 'NR>1 {print $1}' $L2_DIR/L2R1_mapping/L2R1_L1R1.hifislida3.o | sort --parallel=$N_THREADS | uniq | wc -l)

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
m15=$(grep -w "Number of input reads" $L2_DIR/L2R2_mapping/genome/L2R2_genome.Log.final.out |cut -d "|" -f2 | sed 's/\t//g')

a=$(echo "scale=4 ; $m15 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m16=$a"%"


########## GENOME

##### Number of HiFi-Slide L2R2 uniquely mapped to genome
# grep -w "Uniquely mapped reads number" $L2_DIR/L2R2_mapping/genome/L2R2_genome.Log.final.out | cut -d "|" -f2 | sed 's/\t//g'

##### Number of HiFi-Slide L2R2 uniquely mapped to genome and to annotated genes
m17=$(cut -f4 $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m17 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m18=$a"%"

##### Number of HiFi-Slide read pairs genome mapped and spatially resolved
m19=$(awk 'FNR==NR {a[$1]; next} FNR> 0 && $4 in a' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.o $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed | cut -f4 | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m19 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m20=$a"%"


##### Number of HiFi reads spatially resolved, under ROI and aligned to genome
m21=$(awk -F '\t' 'NR>1 {print $1}' $L2_DIR/L2R1_L2R2_integrate/HiFi_L2R2_genome_spatial.txt | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m21 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m22=$a"%"

# Number of genes 
m23=$(awk -F '\t' 'NR>1 {print $9}' $L2_DIR/L2R1_L2R2_integrate/HiFi_L2R2_genome_spatial.txt | sort --parallel=$N_THREADS | uniq | wc -l)

##### Average number of HiFi-Slide read pairs genome mapped and spatially resolved per tile
n_tiles=$(wc -l $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o | cut -d " " -f 1)
m24=$(echo "scale=4 ; $m19 / $n_tiles" | bc | awk '{printf("%.2f",$1)}')

##### Average number of HiFi-Slide read pairs genome mapped and spatially resolved per tile under ROI
n_tiles_under_ROI=$(wc -l $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt | cut -d " " -f 1)
m25=$(echo "scale=4 ; $m21 / $n_tiles_under_ROI" | bc | awk '{printf("%.2f",$1)}')


#### Labels
M1="Total number of barcodes"
M2="Number of deduplicated barcodes"
M3="Percentage of deduplicated barcodes"
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
for k in $(seq 1 25); do
Mk=M${k}
mk=m${k}
echo -e ${!Mk}'\t'${!mk}>> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
done


########## TRANSCRIPTOME
m27=$m19 # to store number of read pairs genome and transcriptome mapped and spatially resolved
m29=$m21 # to store number of read pairs genome and transcriptome mapped and undr ROI

for my_transcript in tRNA piRNA miRNA circRNA; do

size=$(stat -c %s $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt")

if [ $size != 0 ]; then

##### Number of uniquely mapped reads to the transcriptome
m26a=$(grep -w "aligned exactly 1 time" $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_mapped.log" | cut -d " " -f5)

a=$(echo "scale=4 ; $m26a / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m26b=$a"%"

##### Number of HiFi-Slide read pairs transcriptome mapped and spatially resolved
m26c=$(awk 'FNR==NR {a[$1]; next} FNR> 0 && $1 in a' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.o $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt" | cut -f1 | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m26c / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m26d=$a"%"

m27=$(($m27 + $m26c)) # update value

##### Number of HiFi reads spatially resolved, under ROI and aligned to transcriptome
m26e=$(awk -F '\t' 'NR>1 {print $1}' $L2_DIR/L2R1_L2R2_integrate/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m26c / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m26f=$a"%"

m29=$(($m29 + $m26e)) # update value

# Number of transcripts 
m26g=$(awk -F '\t' 'NR>1 {print $6}' $L2_DIR/L2R1_L2R2_integrate/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" | sort --parallel=$N_THREADS | uniq | wc -l)

##### Average number of HiFi-Slide read pairs transcriptome mapped and spatially resolved per tile
n_tiles=$(wc -l $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o | cut -d " " -f 1)
m26h=$(echo "scale=4 ; $m26c / $n_tiles" | bc | awk '{printf("%.2f",$1)}')

##### Average number of HiFi-Slide read pairs transcriptome mapped and spatially resolved per tile under ROI
n_tiles_under_ROI=$(wc -l $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt | cut -d " " -f 1)
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




















