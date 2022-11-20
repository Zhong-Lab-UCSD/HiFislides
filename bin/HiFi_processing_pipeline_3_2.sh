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


### FASTP filtering
fastp \
-i $L2_DIR/L2R2_preprocessing/L2R2.pear_filter.fastq \
-o $L2_DIR/L2R2_preprocessing/L2R2.$RUNNING_LABEL.fastq \
-h $L2_DIR/L2R2_preprocessing/L2R2.$RUNNING_LABEL.log.html \
-j $L2_DIR/L2R2_preprocessing/L2R2.$RUNNING_LABEL.log.json \
--trim_poly_g \
--trim_poly_x \
--disable_quality_filtering \
--thread 16

echo "[$(date '+%m-%d-%y %H:%M:%S')] Filtering complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Align HiFi R2 reads to genome/genes in order to obtain gene annotation for HiFi read pairs.
echo "[$(date '+%m-%d-%y %H:%M:%S')] Align HiFi-Slide reads R2 (L2R2) to the genome using STAR..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_fastp_filter
L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.genome_fastp_filter.fastq

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

L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping/genome_nofilter
L2R2_FILTER_FASTQ=$L2R2_FASTQ

mkdir -p $L2R2_GENOME_DIR
STAR \
--genomeDir $STAR_INDEX \
--readFilesIn $L2R2_FILTER_FASTQ \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx \
--outSAMattributes All \
--outFileNamePrefix $L2R2_GENOME_DIR/L2R2_genome. \
--sjdbGTFfile $annotation_gtf_file \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--runThreadN $N_THREADS


### Select uniquely mapped reads
# samtools view -@ $N_THREADS -b -h -q 255 \
# -o $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
# $L2R2_GENOME_DIR/L2R2_genome.Aligned.sortedByCoord.out.bam

samtools view -@ $N_THREADS -b -h -q 255 \
-o $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
$L2R2_GENOME_DIR/L2R2_genome.Aligned.out.bam

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

### Reads not overlapping annotated genes 
bedtools intersect \
-a $L2R2_GENOME_DIR/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
-b $annotation_gtf_file_genes \
-v -bed | cut -f 4 | awk '{print $0, "\tNA", "\tNA", "\tNA"}' > $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed

### Concatenate and sort
cat $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_1.bed $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed | sort -k 1 --parallel=$N_THREADS -S 20G > $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed

# rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_1.bed
# rm $L2R2_GENOME_DIR/HiFi_L2R2_genome_temp_2.bed

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
L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping/$RUNNING_LABEL
mkdir -p $L2R1_MAPPING_DIR

# Calculate the minimum base match (-k) as 50% of the length of the spatial barcodes (TBD)
temp=$(less $L1R1_DEDUP | head -2 | sed -n '2p')
min_base_match=$(echo "scale=4 ; ${#temp} * 0.5" | bc | awk '{printf("%.0f",$1)}')

bwa mem \
-a \
-k $min_base_match \
-t $N_THREADS \
$L1R1_FASTQ_BWA_INDEX $L2R1_FASTQ > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam 2>$L2R1_MAPPING_DIR/L2R1_L1R1_dedup_$min_base_match.log

### Remove header (not useful and it only occupies storage)
grep -v '^@' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup_k$min_base_match.sam

rm $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam

echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment done."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Alignment done." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

### Select flowcell surface if flowcell type is MiniSeq
if [ "$flowcell_type" == "MiniSeq" ]; then
echo "[$(date '+%m-%d-%y %H:%M:%S')] Selecting the correct flowcell surface..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

# Select 0 and 256 flags
awk -F"\t" '$2 == "0" || $2 == "256" { print $0 }' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup_k$min_base_match.sam > $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.temp.sam

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

L2R1_L1R1_SAM=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup_k$min_base_match.sam

fi


### Filter SAM file to select only HiFi-Slide reads mapped to genome/transcriptome (samf: "SAM filter" custom format)
L2R1_L1R1_SAM=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/IGM/HiFi_placenta_1/L2R1_mapping/L2R1_L1R1_dedup_k19.sam

L2R1_L1R1_SAM_FILTER=$L2R1_MAPPING_DIR/L2R1_L1R1_dedup.filter.sam

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed $L2R1_L1R1_SAM > $L2R1_L1R1_SAM_FILTER

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

echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location (hifislida3.pl) complete."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Match HiFi-Slide R1 reads under ROI with their spatial location (hifislida3.pl) complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide L2R1 complete."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Processing HiFi-Slide L2R1 complete." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log




########## Integrate spatial coordinates and gene expression information
echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information..." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log

L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate/$RUNNING_LABEL
mkdir -p $L2R1_L2R2_INTEGRATE_DIR

### Genome
HiFi_L2R2_genome=$L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed

join -1 1 -2 1 -t $'\t' $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort.o $HiFi_L2R2_genome | cut -f 2,3,4,5,6,7,8 > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt

# rm $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.txt

# Calculate the number of spots in each tile
awk -v OFS='\t' '{print $1, $1"_"$2"_"$3}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt | awk '!seen[$2]++' | awk '{A[$1]++}END{for(i in A)print i,A[i]}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_spot_number_table.txt

# Calculate the number of genes in each tile
awk -v OFS='\t' '{print $1, $1"_"$5}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt | awk '!seen[$2]++' | awk '{A[$1]++}END{for(i in A)print i,A[i]}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_gene_number_table.txt

# Calculate the total gene expression in each tile
awk -v OFS='\t' '{print $1, $7}' $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt | awk -F"\t" '{array[$1]+=$2} END { for (i in array) {print i"\t" array[i]}}' | awk -v OFS='\t' '{print $1, $2, $2/10000}' > $L2R1_L2R2_INTEGRATE_DIR/tile_gene_expression_table.txt


echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Integrate spatial coordinates and gene expression information finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


########## Select spots in ROI after visual inspection of the image
ROI_TILES=$L2_DIR/ROI_tiles_1.txt
ROI_label="ROI_1"

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1"_"$2"_"$3, $1, $2, $3, $5, $6, $7}' | awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]+=$7} END { for (i in array) {print i"\t" array[i]}}' > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.$ROI_label.txt

ROI_TILES=$L2_DIR/ROI_tiles_2.txt
ROI_label="ROI_2"

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1"_"$2"_"$3, $1, $2, $3, $5, $6, $7}' | awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]+=$7} END { for (i in array) {print i"\t" array[i]}}' > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.$ROI_label.txt

ROI_TILES=$L2_DIR/ROI_tiles_3.txt
ROI_label="ROI_3"

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1"_"$2"_"$3, $1, $2, $3, $5, $6, $7}' | awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]+=$7} END { for (i in array) {print i"\t" array[i]}}' > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.$ROI_label.txt

ROI_TILES=$L2_DIR/tile_max_spot.txt
ROI_label="tile_max_spot"

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1"_"$2"_"$3, $1, $2, $3, $5, $6, $7}' | awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]+=$7} END { for (i in array) {print i"\t" array[i]}}' > $L2R1_L2R2_INTEGRATE_DIR/HiFi_L2R2_genome_spatial.final.$ROI_label.txt



####################### QC metrics
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

##### Number of HiFi-Slide L2R2 passing length filtering (performed automatically by fastp)
# a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq | cut -d " " -f 1)
# echo $(($a / 4))
# To make it faster, we can just count the number of input reads in the log of the STAR aligner
m4=$(grep -w "Number of input reads" $L2R2_GENOME_DIR/L2R2_genome.Log.final.out |cut -d "|" -f2 | sed 's/\t//g')

a=$(echo "scale=4 ; $m4 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m5=$a"%"

##### Number of HiFi-Slide L2R2 uniquely mapped to genome and to annotated genes
m6=$(cut -f4 $L2R2_GENOME_DIR/HiFi_L2R2_genome_ALL.sort.bed | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m6 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m7=$a"%"

##### Number of HiFi read pairs mapped to the genome and spatially resolved
m8=$(wc -l $L2R1_MAPPING_DIR/L2R1_L1R1.hifislida3.sort.o | cut -d " " -f 1)

a=$(echo "scale=4 ; $m8 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m9=$a"%"


##### Print to file
M1="Number of HiFi-Slide read pairs"
M2="Number of HiFi-Slide read pairs passing PEAR filtering"
M3="Percentage of HiFi-Slide read pairs passing PEAR filtering"
M4="Number of HiFi-Slide read pairs passing PEAR and FASTP filtering"
M5="Percentage of HiFi-Slide read pairs passing PEAR and FASTP filtering"
M6="Number of HiFi-Slide read pairs genome mapped"
M7="Percentage of HiFi-Slide read pairs genome mapped"
M8="Number of HiFi-Slide read pairs genome mapped and spatially resolved"
M9="Percentage of HiFi-Slide read pairs genome mapped and spatially resolved"


rm $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
touch $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".QC_metrics.txt"
for k in $(seq 1 7); do
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




















