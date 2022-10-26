################## INPUT PARAMETERS
BIN_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/bin

N_THREADS=32
BWA_MEMORY=80000 # memory (in Megabytes) to be used for bwa index. It does not seem that useful.

# Flowcell and surface identifiers
flowcell_type="NextSeq" # one of: MiniSeq, NextSeq
flowcell="AAAL33WM5"

if [ "$flowcell_type" == "NextSeq" ]; then
surface=$flowcell:1:1
elif [ "$flowcell_type" == "MiniSeq" ]; then
surface=$flowcell:1:
fi

# Directories of the processed data
L1_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/$flowcell # spatial barcodes
mkdir -p $L1_DIR

# Directories of the raw fastq files for each library. The full path is used here.
L1_FASTQ_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/MiniSeq/test_sample/lib1/fastq
L1_FASTQ_BASENAME=*_L001_R1_001.fastq.gz


################## PROCESSING
touch $L1_DIR/$flowcell.log

echo "Processing of "$flowcell
echo "Processing of "$flowcell >> $L1_DIR/$flowcell.log
echo "------------------------------" >> $L1_DIR/$flowcell.log


#################### LIBRARY 1 (spatial barcodes)
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L1R1..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L1R1..." >> $L1_DIR/$flowcell.log
echo "------------------------------" >> $L1_DIR/$flowcell.log

### Deduplication of raw reads from the recycled flow cell to extract unique raw reads as spatial barcodes
# g++ surfdedup.cpp -o surfdedup -lz
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start deduplication of L1R1 reads..."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start deduplication of L1R1 reads..." >> $L1_DIR/$flowcell.log

$BIN_DIR/surfdedup \
$surface \
$L1_FASTQ_DIR/$L1_FASTQ_BASENAME > $L1_DIR/L1R1_dedup.fasta 2>$L1_DIR/L1R1_dup.txt

L1R1_DEDUP=$L1_DIR/L1R1_dedup.fasta
L1R1_DUP=$L1_DIR/L1R1_dup.txt

echo "[$(date '+%m-%d-%y %H:%M:%S')] Deduplication of L1R1 reads complete."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Deduplication of L1R1 reads complete." >> $L1_DIR/$flowcell.log

# dummy example to make running faster
# $BIN_DIR/surfdedup $surface $L1_FASTQ_DIR/MT080_S1_L001_R1_001.fastq.gz > $L1_DIR/L1R1_dedup.fasta 2>$L1_DIR/L1R1_dup.txt

### Align HiFi R1 reads (L2R1) to spatial barcodes (L1R1) in order to obtain spatial coordinates for HiFi read pairs.

# Create index files for L1R1
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start creating BWA index for spatial barcodes (L1R1)..." >> $L1_DIR/$flowcell.log
BWA_BLOCK_SIZE=$(($BWA_MEMORY * 1000000 / 8)) # currently not used
mkdir -p $L1_DIR/bwa_index_L1R1

bwa index \
-p $L1_DIR/bwa_index_L1R1/L1R1_dedup \
$L1R1_DEDUP

echo "[$(date '+%m-%d-%y %H:%M:%S')] BWA index creation complete." >> $L1_DIR/$flowcell.log

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Finished processing HiFi-Slide L1R1."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Finished processing HiFi-Slide L1R1." >> $L1_DIR/$flowcell.log
echo "------------------------------" >> $L1_DIR/$flowcell.log


####################### QC metrics
echo "------------------------------" >> $L1_DIR/$flowcell.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..." >> $L1_DIR/$flowcell.log

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
a=$(wc -l $L1R1_DEDUP | cut -d " " -f 1)
m2=$(($a / 2))

a=$(echo "scale=4 ; $m2 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m3=$a"%"

#### Labels
M1="Total number of barcodes"
M2="Number of deduplicated barcodes"
M3="Percentage of deduplicated barcodes"

rm $L1_DIR/$flowcell".QC_metrics.txt"
touch $L1_DIR/$flowcell".QC_metrics.txt"
for k in $(seq 1 3); do
Mk=M${k}
mk=m${k}
echo -e ${!Mk}'\t'${!mk}>> $L1_DIR/$flowcell".QC_metrics.txt"
done

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] QC metrics calculation finished." 
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] QC metrics calculation finished." >> $L1_DIR/$flowcell.log






