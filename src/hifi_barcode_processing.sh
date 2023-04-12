#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    ---------------------------------------------------------------------------
    Usage: $PROGNAME [-b <bin_dir>]
                     [-f <flowcell_id>] [-l <flowcell_lane>] [-s <flowcell_surface>]
                     [-d <spatial_barcode_dir>] [-N <name_fastq.gz>] [-T <tiles.txt>]
                     [-t <threads>] [-o <output_dir>]
    
    Dependency: bwa
    This is a script to pre-process spatial barcodes.
    -b : Directory of the scripts.
    -f : Flowcell ID.
    -l : Flowcell lane.
    -s : Flowcell surface (1 or 2).
    -d : Directory of the barcode fastq files.
    -N : Name of the fastq file (also *R1_001.fastq.gz is accepted).
    -T : Tiles to be used, one per line (if empty, all the tiles in the flowcell are used).
    -t : Max CPU threads for parallelized processing, at least 4 (default 8).
    -o : Output directory.
    -h : Show usage help
    ---------------------------------------------------------------------------
EOF
    exit 1
}

parameter_error() {
    echo "!!!!!!!!!!!!!!!! Pipeline exit with parameter ERROR at $(date '+%Y-%m-%d %H:%M:%S %Z') !!!!!!!!!!!!!!!!"
    usage
}

while getopts :b:f:l:s:d:N:T:t:o:h opt; do
    case $opt in
        b) BIN_DIR=${OPTARG};;
        f) flowcell=${OPTARG};;
        l) flowcell_lane=${OPTARG};;
        s) flowcell_surface=${OPTARG};;
        d) L1_FASTQ_DIR=${OPTARG};;
        N) L1_FASTQ_FILENAME=${OPTARG};;
        T) tiles=${OPTARG};;
        t) N_THREADS=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        h) usage;;
    esac
done

# check parameters
[ -z "$BIN_DIR" ] &&  echo "Error!! Please provide path to scripts with -b" && parameter_error

[  -z "$flowcell" ] && echo "Error!! Please provide the flowcell ID with -f" && parameter_error

[  -z "$flowcell_lane" ] && echo "Error!! Please provide the flowcell lane with -l" && parameter_error

[  -z "$flowcell_surface" ] && echo "Error!! Please provide the flowcell surface (1 or 2) with -s" && parameter_error

[ ! -d "$L1_FASTQ_DIR" ] && echo "Error!! Directory of spatial barcode fastq files does not exist: "$L1_FASTQ_DIR && parameter_error

[ -z "$L1_FASTQ_FILENAME" ] && echo "Error!! Filename of the fastq file (you can also use for example *R1_001.fastq.gz if multiple there are multiple files) with -N" && parameter_error

[  -z "$N_THREADS" ] && echo "Use default thread number 8'." && N_THREADS=8
if ! [[ "$N_THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && parameter_error 
fi


################## PROCESSING
FLOWCELL_FULL=$flowcell"_"$flowcell_lane"_"$flowcell_surface

# Directories of the processed data
L1_DIR=$OUT_DIR/$FLOWCELL_FULL
mkdir -p $L1_DIR

> $L1_DIR/$FLOWCELL_FULL.log

echo "Processing of "$FLOWCELL_FULL
echo "Processing of "$FLOWCELL_FULL >> $L1_DIR/$FLOWCELL_FULL.log
echo "------------------------------" >> $L1_DIR/$FLOWCELL_FULL.log


#################### LIBRARY 1 (spatial barcodes)
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L1R1..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start processing HiFi-Slide L1R1..." >> $L1_DIR/$FLOWCELL_FULL.log
echo "------------------------------" >> $L1_DIR/$FLOWCELL_FULL.log

### Deduplication of raw reads from the recycled flow cell to extract unique raw reads as spatial barcodes
# g++ surfdedup.cpp -o surfdedup -lz
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start deduplication of L1R1 reads..."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start deduplication of L1R1 reads..." >> $L1_DIR/$FLOWCELL_FULL.log

if [ -z "$tiles" ]; then
$BIN_DIR/surfdedup \
$flowcell":"$flowcell_lane":"$flowcell_surface \
$L1_FASTQ_DIR/$L1_FASTQ_FILENAME > $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta 2>$L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt
else
$BIN_DIR/surfdedup \
$flowcell":"$flowcell_lane":"$flowcell_surface \
$tiles \
$L1_FASTQ_DIR/$L1_FASTQ_FILENAME > $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta 2>$L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt
fi

echo "[$(date '+%m-%d-%y %H:%M:%S')] Deduplication of L1R1 reads complete."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Deduplication of L1R1 reads complete." >> $L1_DIR/$FLOWCELL_FULL.log

### Align HiFi R1 reads (L2R1) to spatial barcodes (L1R1) in order to obtain spatial coordinates for HiFi read pairs.

# Create index files for L1R1
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start creating BWA index for spatial barcodes (L1R1)..." >> $L1_DIR/$FLOWCELL_FULL.log

mkdir -p $L1_DIR/bwa_index_L1R1

bwa index \
-p $L1_DIR/bwa_index_L1R1/$FLOWCELL_FULL.L1R1_dedup \
$L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta

echo "[$(date '+%m-%d-%y %H:%M:%S')] BWA index creation complete." >> $L1_DIR/$FLOWCELL_FULL.log

### Count the number of barcodes in each tile
echo "[$(date '+%m-%d-%y %H:%M:%S')] Count the number of barcodes in each tile..."

grep ">" $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta | cut -f5 -d ":" | awk -v OFS='\t' '{A[$1]++}END{for(i in A)print i,A[i]}' > $L1_DIR/tile_barcode_number_temp1.txt
awk -v OFS='\t' '{print $2}' $L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt | cut -f5 -d ":" | awk -v OFS='\t' '{A[$1]++}END{for(i in A)print i,A[i]}' > $L1_DIR/tile_barcode_number_temp2.txt

cat $L1_DIR/tile_barcode_number_temp1.txt $L1_DIR/tile_barcode_number_temp2.txt | awk -v OFS='\t' '{a[$1]+=$2}END{for(i in a) print i,a[i]}' > $L1_DIR/tile_barcode_number.txt

rm $L1_DIR/tile_barcode_number_temp1.txt
rm $L1_DIR/tile_barcode_number_temp2.txt

# Uniquely located barcodes
grep "_1" $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta | cut -f5 -d ":" | awk -v OFS='\t' '{A[$1]++}END{for(i in A)print i,A[i]}' > $L1_DIR/tile_barcode_number_N1.txt

echo "[$(date '+%m-%d-%y %H:%M:%S')] Count the number of barcodes in each tile complete."


echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Finished processing HiFi-Slide L1R1."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Finished processing HiFi-Slide L1R1." >> $L1_DIR/$FLOWCELL_FULL.log
echo "------------------------------" >> $L1_DIR/$FLOWCELL_FULL.log


####################### QC metrics
echo "------------------------------" >> $L1_DIR/$FLOWCELL_FULL.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..." >> $L1_DIR/$FLOWCELL_FULL.log

##### Total number of barcodes (L1R1)
> $L1_DIR/$FLOWCELL_FULL.L1R1_stats.txt

for x in $L1_FASTQ_DIR/$L1_FASTQ_FILENAME; do
echo $x
a=$(unpigz -p $N_THREADS -c $x | wc -l)
out=$(($a / 4))
echo -e $x"\t"$out >> $L1_DIR/$FLOWCELL_FULL.L1R1_stats.txt
done

m1=$(awk '{ sum += $2 } END { print sum }' $L1_DIR/$FLOWCELL_FULL.L1R1_stats.txt)

##### Number of deduplicated barcodes
a=$(wc -l $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta | cut -d " " -f 1)
m2=$(($a / 2))

a=$(echo "scale=4 ; $m2 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m3=$a"%"

#### Labels
M1="Total number of barcodes"
M2="Number of deduplicated barcodes"
M3="Percentage of deduplicated barcodes"

> $L1_DIR/$FLOWCELL_FULL".QC_metrics.txt"
for k in $(seq 1 3); do
Mk=M${k}
mk=m${k}
echo -e ${!Mk}'\t'${!mk}>> $L1_DIR/$FLOWCELL_FULL".QC_metrics.txt"
done

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] QC metrics calculation finished." 
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] QC metrics calculation finished." >> $L1_DIR/$FLOWCELL_FULL.log


#######

# /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/src/bin/surfdedup \
# AAANLCHHV:1:1 \
# /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/tiles_top_14x6.txt \
# /mnt/extraids/SDSC_NFS/linpei/hifi/recycled_flowcell/2022_09_06_NS/*R1_001.fastq.gz > /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANLCHHV_1_1/top/AAANLCHHV_1_1.L1R1_dedup.fasta 2>/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANLCHHV_1_1/top/AAANLCHHV_1_1.L1R1_dup.txt

# /mnt/SDSC_NFS/rcalandrelli/HiFi/data/src/bin/surfdedup \
# AAANLCHHV:1:1 \
# /mnt/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/tiles_bottom_14x6.txt \
# /mnt/SDSC_NFS/linpei/hifi/recycled_flowcell/2022_09_06_NS/*R1_001.fastq.gz > /mnt/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANLCHHV_1_1/bottom/AAANLCHHV_1_1.L1R1_dedup.fasta 2>/mnt/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANLCHHV_1_1/bottom/AAANLCHHV_1_1.L1R1_dup.txt


# /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/src/bin/surfdedup \
# AAANLCHHV:2:1 \
# /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/tiles_top_14x6.txt \
# /mnt/extraids/SDSC_NFS/linpei/hifi/recycled_flowcell/2022_09_06_NS/*R1_001.fastq.gz > /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANLCHHV_2_1/top/AAANLCHHV_2_1.L1R1_dedup.fasta 2>/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANLCHHV_2_1/top/AAANLCHHV_2_1.L1R1_dup.txt

# /mnt/SDSC_NFS/rcalandrelli/HiFi/data/src/bin/surfdedup \
# AAANLCHHV:2:1 \
# /mnt/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/tiles_bottom_14x6.txt \
# /mnt/SDSC_NFS/linpei/hifi/recycled_flowcell/2022_09_06_NS/*R1_001.fastq.gz > /mnt/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANLCHHV_2_1/bottom/AAANLCHHV_2_1.L1R1_dedup.fasta 2>/mnt/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANLCHHV_2_1/bottom/AAANLCHHV_2_1.L1R1_dup.txt

### To extract coordinates of all the barcodes (USELESS)
# tile_matrix_expanded=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/tile_matrix_expanded_14x6.txt
# input_file_dup=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAAV7J7HV_1_1/AAAV7J7HV_1_1.L1R1_dup.txt

# input_file_dup=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAAV7J7HV_1_1/temp_dup.txt

# barcode_id=$(awk -v OFS='\t' '{print $2}' $input_file_dup)
# tile_id=$(cut -d':' -f5 $barcode_id)
# col=$(cut -d':' -f6 $barcode_id)
# row=$(cut -d':' -f7 $barcode_id)
# paste <(echo "$barcode_id") <(echo "$tile_id") <(echo "$col") <(echo "$row") --delimiters '\t' | sort -k 2 > temp.txt


# join -1 3 -2 2 -o 2.1,2.2,2.3,2.4,1.4,1.5 <(sort -k 3 $tile_matrix_expanded) <(sort -k 2 temp.txt) -t $'\t' | awk -v OFS='\t' '{print $1, $3+$5, $4+$6}' > AAAV7J7HV_1_1.L1R1_dup.coord.txt
