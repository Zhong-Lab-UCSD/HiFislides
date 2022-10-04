BIN_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/bin

OUT_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data
SAMPLE_NAME=test_sample_new

# Directories of the raw fastq files for each library. The full path is used here.
L1_FASTQ_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/test_sample/lib1/fastq
L2_FASTQ_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/test_sample/lib2/fastq

# Directories of the processed data
L1_DIR=$OUT_DIR/$SAMPLE_NAME/lib1 # spatial barcodes
L2_DIR=$OUT_DIR/$SAMPLE_NAME/lib2 # HiFi library
mkdir -p $L1_DIR
mkdir -p $L2_DIR

# Raw reads of HiFi Slides sequencing
L2R1=$L2_FASTQ_DIR/Undetermined_S0_L001_R1_001.fastq.gz
L2R2=$L2_FASTQ_DIR/Undetermined_S0_L001_R2_001.fastq.gz

# Flowcell and surface identifiers
# Q: should i and j be elevated as input parameters and kept into the filenames as currently we are doing? Discuss further.
flowcell=AAAL33WM5
i=1
j=1
surf=$flowcell:$i:$j

### Annotation file hg38. This file can be downloaded without the need of computing it from the GTF file or bedtools intersect can take GTF as input, only genes can be selected from the GTF file.

# annotation_gtf_file=/dataOS/sysbio/Genomes/Homo_sapiens/Ensembl/GRCH38_hg38/Annotation/Genes/Homo_sapiens.GRCh38.84.chr.gtf
annotation_gtf_file=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gtf

# Select full genes only

# awk -v OFS='\t' '$3=="gene"' $annotation_gtf_file > /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/Homo_sapiens.GRCh38.84.chr.gene.gtf
awk -v OFS='\t' '$3=="gene"' $annotation_gtf_file > /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gene.gtf

# awk -v OFS='\t' '$3=="exon"' $annotation_gtf_file > /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/Homo_sapiens.GRCh38.84.chr.gene.gtf # to select transcripts only?

### Mapping reference indexes
STAR_INDEX=/dataOS/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/STARindex_withSJ
BOWTIE2_INDEX=/mnt/extraids/SDSC_NFS/linpei/genome/HSATR




########## LIBRARY 1 (spatial barcodes)

### Deduplication of raw reads from the recycled flow cell to extract unique raw reads as spatial barcodes

# g++ surfdedup.cpp -o surfdedup -lz
$BIN_DIR/surfdedup $surf $L1_FASTQ_DIR/*_L00$i"_R1_001.fastq.gz" > $L1_DIR/L1R1_dedup_$i"_"$j.fasta 2>$L1_DIR/L1R1_dup_$i"_"$j.txt

# dummy example to make running faster
# $BIN_DIR/surfdedup $surf $L1_FASTQ_DIR/MT080_S1_L001_R1_001.fastq.gz > $L1_DIR/L1R1_dedup_$i"_"$j.fasta 2>$L1_DIR/L1R1_dup_$i"_"$j.txt


### Align HiFi R1 reads (L2R1) to spatial barcodes (L1R1) in order to obtain spatial coordinates for HiFi read pairs.

# Create index files for L2R1
mkdir -p $L1_DIR/bwa_index_L1R1
# cd $L1_DIR/bwa_index_L1R1
bwa index -p $L1_DIR/bwa_index_L1R1/L1R1_dedup_$i"_"$j $L1_DIR/L1R1_dedup_$i"_"$j.fasta

# Alignment
mkdir -p $L2_DIR/L2R1_mapping
bwa mem -a -k 40 -t 32 $L1_DIR/bwa_index_L1R1/L1R1_dedup_$i"_"$j $L2R1 > $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".sam" 2>$L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".log"

### Select HiFi-Slide R1 reads with highest alignment score
$BIN_DIR/hifislida.pl $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".sam" > $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".hifislida.o" 2>$L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".hifislida.e"

### Rank the tiles by number of HiFi-Slide read pairs
$BIN_DIR/hifislida2.pl \
$L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".hifislida.o" \
$L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".sam" > $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".hifislida2.o" 2>$L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".hifislida2.e"

### Select tiles under ROI
$BIN_DIR/select_tiles_in_ROI.r \
-i $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".hifislida2.o" \
-o $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt \
--max_size_ROI 4 \
--min_size_ROI 2 \
--p_value 0.05

### Match HiFi-Slide read pairs with spatial location
$BIN_DIR/hifislida3.pl \
$L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_$i"_"$j".hifislida.o" \
$L2_DIR/L2R1_mapping/ROI_tile_IDs.txt \
$L1_DIR/L1R1_dup_$i"_"$j.txt > $L2_DIR/L2R1_mapping/temp.hifislida3.o

# Add header
echo -e "HiFi_read_id\ttile_id\tcol\trow\tN" | cat - $L2_DIR/L2R1_mapping/temp.hifislida3.o > $L2_DIR/L2R1_mapping/L2R1__L1R1_$i"_"$j".hifislida3.o"

rm $L2_DIR/L2R1_mapping/temp.hifislida3.o



########## LIBRARY 2 (HiFi-Slide read pairs)

### Preprocessing of HiFi R2 reads
mkdir -p $L2_DIR/L2R2_preprocessing

# Find L2R2 reads overlapping L2R1 and filter them out using the software pear
minoverlap=10

pear \
-f $L2R1 \
-r $L2R2 \
-v $minoverlap \
-j 32 \
-o $L2_DIR/L2R2_preprocessing/L2R2_pear

# grep @ $L2_DIR/L2R2_preprocessing/L2R2_pear.unassembled.reverse.fastq > $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.names
# zgrep -A 3 -f $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.names $L2R2 > $L2_DIR/L2R2_preprocessing/L2R2.pear_filter_names.fastq

# Necessary? To be confirmed!
seqtk seq -r $L2_DIR/L2R2_preprocessing/L2R2_pear.unassembled.reverse.fastq > $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq

# Trimming the front 60 bp to remove Illumina adapters (max 16 threads allowed, max 30 bp at the time)
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


### Align HiFi R2 reads to genome/genes in order to obtain gene annotation for HiFi read pairs.

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
--runThreadN 32

### Select uniquely mapped reads
samtools view -@ 32 -b -h -q 255 \
-o $L2_DIR/L2R2_mapping/genome/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
$L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam

samtools view -@ 32 $L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam | wc -l

samtools view -@ 32 -q 255 $L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam | wc -l

samtools view -@ 32 -q 30 $L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam | wc -l


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


### Align HiFi R2 reads to the transcriptome in order to obtain gene annotation for HiFi read pairs.
mkdir -p $L2_DIR/L2R2_mapping/transcriptome

### Creating Bowtie 2 indexes
for t in tRNA piRNA mirbase circbase; do
mkdir -p /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/$t/bowtie2_index

bowtie2-build \
--threads 32 \
/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/$t/*fa* \
/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/$t/bowtie2_index/$t
done

# Dummy fastq for testing purposes
# head -40000 $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq > $L2_DIR/L2R2_mapping/transcriptome/temp_L2R2.trim_front_60.fastq

for t in tRNA piRNA mirbase circbase; do
mkdir -p $L2_DIR/L2R2_mapping/transcriptome/$t

bowtie2 \
-x /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/$t/bowtie2_index/$t \
-U $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq \
-S $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.sam" \
--un $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_unmapped.txt" \
--no-unal --threads 32 --local 2> $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.log"

# Select uniquely mapped reads
samtools view -q 10 \
-o $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.mapq10.sam" \
$L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.sam"

# Extract fields of interest
cut -f 1,3 $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.mapq10.sam" > $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.mapq10.txt"

done


# samtools view -q 30 $L2_DIR/L2R2_mapping/transcriptome/L2R2_transcriptome_mapped.sam | wc -l

# samtools view -q 10 \
# -o $L2_DIR/L2R2_mapping/transcriptome/L2R2_transcriptome_mapped.mapq10.sam \
# $L2_DIR/L2R2_mapping/transcriptome/L2R2_transcriptome_mapped.sam

# cut -f 1,3 $L2_DIR/L2R2_mapping/transcriptome/L2R2_transcriptome_mapped.mapq10.sam > $L2_DIR/L2R2_mapping/transcriptome/L2R2_transcriptome_mapped.mapq10.temp.txt




########## Integrate spatial coordinates and gene expression information
mkdir -p $L2_DIR/L2R1_L2R2_integrate

HiFi_L2R1_spatial=L2R1__L1R1_$i"_"$j".hifislida3.sort.o"
#sort -k 1 $L2_DIR/L2R1_mapping/hifislida3.o > $HiFi_L2R1_spatial
cat L2R1__L1R1_$i"_"$j".hifislida3.o" | sort -k 1 --parallel=32 -S 20G > $HiFi_L2R1_spatial

### Genome
HiFi_L2R2_genome=$L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.sort.bed
#sort -k 4 $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed > $HiFi_L2R2_genome
cat $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed | sort -k 4 --parallel=32 -S 20G > $HiFi_L2R2_genome

join -1 1 -2 4 -t $'\t' $HiFi_L2R1_spatial $HiFi_L2R2_genome > $L2_DIR/L2R1_L2R2_integrate/temp_HiFi_L2R2_genome_spatial.txt

# Add header
echo -e "HiFi_read_id\ttile_id\tcol\trow\tN\tHiFi_read_chr\tHiFi_read_start\tHiFi_read_end\tgene_id\tgene_name\tgene_type" | cat - $L2_DIR/L2R1_L2R2_integrate/temp_HiFi_L2R2_genome_spatial.txt > $L2_DIR/L2R1_L2R2_integrate/HiFi_L2R2_genome_spatial.txt

rm $L2_DIR/L2R1_L2R2_integrate/temp_HiFi_L2R2_genome_spatial.txt

### Transcriptome

for t in tRNA piRNA mirbase circbase; do

size=$(stat -c %s $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.mapq10.txt")

if [ $size != 0 ]; then

mkdir -p $L2_DIR/L2R1_L2R2_integrate/$t

cat $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.mapq10.txt" | sort -k 1 --parallel=32 -S 20G > $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.mapq10.sort.txt"

join -1 1 -2 1 -t $'\t' $HiFi_L2R1_spatial $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.mapq10.sort.txt" > $L2_DIR/L2R1_L2R2_integrate/$t/temp_HiFi_L2R2_$t"_spatial.txt"

# Add header
echo -e "HiFi_read_id\ttile_id\tcol\trow\tN\ttranscript_id" | cat - $L2_DIR/L2R1_L2R2_integrate/$t/temp_HiFi_L2R2_$t"_spatial.txt" > $L2_DIR/L2R1_L2R2_integrate/$t/HiFi_L2R2_$t"_spatial.txt"

rm $L2_DIR/L2R1_L2R2_integrate/$t/temp_HiFi_L2R2_$t"_spatial.txt"

fi

done








