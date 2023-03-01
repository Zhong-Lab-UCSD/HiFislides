#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    ---------------------------------------------------------------------------
    Usage: $PROGNAME [-b <bin_dir>]
                     [-f <flowcell_id>] [-l <flowcell_lane>] [-s <flowcell_surface>]
                     [-d <spatial_barcode_dir>] [-N <suffix_fastq.gz>] [-T <tiles.txt>]
                     [-t <threads>] [-o <output_dir>]
    
    Dependency: bwa
    This is a script to pre-process spatial barcodes.
    -b : Directory of the scripts.
    -f : Flowcell ID.
    -l : Flowcell lane.
    -s : Flowcell surface (1 or 2).
    -d : Directory of the barcode fastq files.
    -N : Suffix of the fastq file (example R1_001.fastq.gz).
    -T : Tiles to be used (if empty, all the tiles in the flowcell are used).
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
        N) L1_FASTQ_SUFFIX=${OPTARG};;
        T) tiles=${OPTARG};;
        t) N_THREADS=${OPTARG};;
        o) L1_DIR=${OPTARG};;
        h) usage;;
    esac
done

# check parameters
[ -z "$BIN_DIR" ] &&  echo "Error!! Please provide path to scripts with -b" && parameter_error

[  -z "$flowcell" ] && echo "Error!! Please provide the flowcell ID with -f" && parameter_error

[  -z "$flowcell_lane" ] && echo "Error!! Please provide the flowcell lane with -l" && parameter_error

[  -z "$flowcell_surface" ] && echo "Error!! Please provide the flowcell surface (1 or 2) with -s" && parameter_error

[ ! -d "$L1_FASTQ_DIR" ] && echo "Error!! Directory of spatial barcode fastq files does not exist: "$L1_FASTQ_DIR && parameter_error

[ -z "$L1_FASTQ_SUFFIX" ] && echo "Error!! Suffix of the fastq file (example R1_001.fastq.gz) with -N" && parameter_error

[  -z "$N_THREADS" ] && echo "Use default thread number 8'." && N_THREADS=8
if ! [[ "$N_THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && parameter_error 
fi


################## PROCESSING
FLOWCELL_FULL=$flowcell"_"$flowcell_lane"_"$flowcell_surface

# Directories of the processed data
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
$L1_FASTQ_DIR/$L1_FASTQ_SUFFIX > $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta 2>$L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt
else
$BIN_DIR/surfdedup \
$flowcell":"$flowcell_lane":"$flowcell_surface \
$tiles \
$L1_FASTQ_DIR/$L1_FASTQ_SUFFIX > $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta 2>$L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt
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

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Finished processing HiFi-Slide L1R1."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Finished processing HiFi-Slide L1R1." >> $L1_DIR/$FLOWCELL_FULL.log
echo "------------------------------" >> $L1_DIR/$FLOWCELL_FULL.log


####################### QC metrics
echo "------------------------------" >> $L1_DIR/$FLOWCELL_FULL.log
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..."
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start QC metrics calculation..." >> $L1_DIR/$FLOWCELL_FULL.log

##### Total number of barcodes (L1R1)
> $L1_DIR/$FLOWCELL_FULL.L1R1_stats.txt

for x in $L1_FASTQ_DIR/$L1_FASTQ_SUFFIX; do
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

rm $L1_DIR/$FLOWCELL_FULL".QC_metrics.txt"
> $L1_DIR/$FLOWCELL_FULL".QC_metrics.txt"
for k in $(seq 1 3); do
Mk=M${k}
mk=m${k}
echo -e ${!Mk}'\t'${!mk}>> $L1_DIR/$FLOWCELL_FULL".QC_metrics.txt"
done

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] QC metrics calculation finished." 
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] QC metrics calculation finished." >> $L1_DIR/$FLOWCELL_FULL.log


######## Test flowcell splitting
n_row_flowcell=14
n_col_flowcell=6

if [ $(($n_row_flowcell%2)) -eq 0 ]; then
end_row_top=$(($n_row_flowcell / 2))
start_row_bottom=$(($n_row_flowcell / 2 + 1))
elif [ $(($n_row_flowcell%2)) -eq 1 ]; then
end_row_top=$(($n_row_flowcell / 2))
start_row_bottom=$(($n_row_flowcell / 2 + 2))
fi

# Tiles in the top portion of the flowcell
> $L1_DIR/tiles_top.txt
for col in $(seq 1100 100 $((1000 + $n_col_flowcell * 100))); do
for row in $(seq 1 $end_row_top); do
echo $(($col+$row)) >> $L1_DIR/tiles_top.txt
done
done

cut -f5 -d ":" $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta | grep -A 1 -n -f $L1_DIR/tiles_top.txt | sed "s/\-/\:/" | cut -f1 -d: | awk 'FNR==NR{h[$1]; c++; next} FNR in h{print; if (!--c) exit}' - $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta > $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.top.fasta

mkdir -p $L1_DIR/bwa_index_L1R1_top

bwa index \
-p $L1_DIR/bwa_index_L1R1_top/$FLOWCELL_FULL.L1R1_dedup.top \
$L1_DIR/$FLOWCELL_FULL.L1R1_dedup.top.fasta

# Tiles in the bottom portion of the flowcell
> $L1_DIR/tiles_bottom.txt
for col in $(seq 1100 100 $((1000 + $n_col_flowcell * 100))); do
for row in $(seq $start_row_bottom $n_row_flowcell); do
echo $(($col+$row)) >> $L1_DIR/tiles_bottom.txt
done
done

cut -f5 -d ":" $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta | grep -A 1 -n -f $L1_DIR/tiles_bottom.txt | sed "s/\-/\:/" | cut -f1 -d: | awk 'FNR==NR{h[$1]; c++; next} FNR in h{print; if (!--c) exit}' - $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta > $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.bottom.fasta

mkdir -p $L1_DIR/bwa_index_L1R1_bottom

bwa index \
-p $L1_DIR/bwa_index_L1R1_bottom/$FLOWCELL_FULL.L1R1_dedup.bottom \
$L1_DIR/$FLOWCELL_FULL.L1R1_dedup.bottom.fasta

#######
cut -f1 $L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt | cut -f5 -d ":" | grep -n -f $L1_DIR/tiles_bottom.txt | cut -f1 -d: > $L1_DIR/dup_bottom1.txt

cut -f2 $L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt | cut -f5 -d ":" | grep -n -f $L1_DIR/tiles_bottom.txt | cut -f1 -d: > $L1_DIR/dup_bottom2.txt

comm -12 $L1_DIR/dup_bottom1.txt $L1_DIR/dup_bottom2.txt > $L1_DIR/dup_bottom.txt

awk 'FNR==NR{h[$1]; c++; next} FNR in h{print; if (!--c) exit}' $L1_DIR/dup_bottom.txt $L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt > $L1_DIR/$FLOWCELL_FULL.L1R1_dup.bottom.txt


##### Count the number of barcodes in each tile
grep ">" $L1_DIR/$FLOWCELL_FULL.L1R1_dedup.fasta | cut -f5 -d ":" | awk '{A[$1]++}END{for(i in A)print i,A[i]}' > $L1_DIR/tile_barcode_number_1.txt
awk -v OFS='\t' '{print $2}' $L1_DIR/$FLOWCELL_FULL.L1R1_dup.txt | cut -f5 -d ":" | awk '{A[$1]++}END{for(i in A)print i,A[i]}' > $L1_DIR/tile_barcode_number_2.txt



