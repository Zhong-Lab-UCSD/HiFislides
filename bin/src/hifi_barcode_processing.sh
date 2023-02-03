#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    ---------------------------------------------------------------------------
    Usage: $PROGNAME [-b <bin_dir>]
                     [-F <flowcell_id>] [-f <flowcell_type>] 
                     [-d <spatial_barcode_dir>] [-N <suffix_fastq.gz>]
                     [-t <threads>] [-o <output_dir>] 
    
    Dependency: bwa
    This is an all-in-one wrapper script of HiFi Pipeline.
    -b : Directory of the scripts.
    -F : Flowcell ID.
    -f : Flowcell type (NextSeq, MiniSeq).
    -d : Directory of the barcode fastq files.
    -N : Suffix of the fastq file (example R1_001.fastq.gz).
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

while getopts :b:F:f:d:N:t:o:h opt; do
    case $opt in
        b) BIN_DIR=${OPTARG};;
        F) flowcell=${OPTARG};;
        f) flowcell_type=${OPTARG};;
        d) L1_FASTQ_DIR=${OPTARG};;
        N) L1_FASTQ_SUFFIX=${OPTARG};;
        t) N_THREADS=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        h) usage;;
    esac
done

# check parameters
[ -z "$BIN_DIR" ] &&  echo "Error!! Please provide path to scripts with -b" && parameter_error

[  -z "$flowcell" ] && echo "Error!! Please provide the flowcell ID with -F" && parameter_error

[  -z "$flowcell_type" ] && echo "Error!! Please provide the flowcell type (NextSeq, MiniSeq) with -f" && parameter_error

[ ! -d "$L1_FASTQ_DIR" ] && echo "Error!! Directory of spatial barcode fastq files does not exist: "$L1_FASTQ_DIR && parameter_error

[ -z "$L1_FASTQ_SUFFIX" ] && echo "Error!! Suffix of the fastq file (example R1_001.fastq.gz) with -N" && parameter_error

[  -z "$N_THREADS" ] && echo "Use default thread number 8'." && N_THREADS=8
if ! [[ "$N_THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && parameter_error 
fi

[ ! -d "$OUT_DIR" ] && echo "Error!! Output directory not exist: "$OUT_DIR && parameter_error


################## PROCESSING

# Directories of the processed data
L1_DIR=$OUT_DIR/$flowcell # spatial barcodes
mkdir -p $L1_DIR

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
$L1_FASTQ_DIR/$L1_FASTQ_SUFFIX > $L1_DIR/$flowcell.L1R1_dedup.fasta 2>$L1_DIR/$flowcell.L1R1_dup.txt

L1R1_DEDUP=$L1_DIR/$flowcell.L1R1_dedup.fasta
L1R1_DUP=$L1_DIR/$flowcell.L1R1_dup.txt

echo "[$(date '+%m-%d-%y %H:%M:%S')] Deduplication of L1R1 reads complete."
echo "[$(date '+%m-%d-%y %H:%M:%S')] Deduplication of L1R1 reads complete." >> $L1_DIR/$flowcell.log

### Align HiFi R1 reads (L2R1) to spatial barcodes (L1R1) in order to obtain spatial coordinates for HiFi read pairs.

# Create index files for L1R1
echo "[$(date '+%m-%d-%y %H:%M:%S')] Start creating BWA index for spatial barcodes (L1R1)..." >> $L1_DIR/$flowcell.log
BWA_BLOCK_SIZE=$(($BWA_MEMORY * 1000000 / 8)) # currently not used
mkdir -p $L1_DIR/bwa_index_L1R1

bwa index \
-p $L1_DIR/bwa_index_L1R1/$flowcell.L1R1_dedup \
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
rm $L1_DIR/$flowcell.L1R1_stats.txt
touch $L1_DIR/$flowcell.L1R1_stats.txt

for x in $L1_FASTQ_DIR/$L1_FASTQ_SUFFIX; do
echo $x
a=$(unpigz -p $N_THREADS -c $x | wc -l)
out=$(($a / 4))
echo -e $x"\t"$out >> $L1_DIR/$flowcell.L1R1_stats.txt
done

m1=$(awk '{ sum += $2 } END { print sum }' $L1_DIR/$flowcell.L1R1_stats.txt)

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




















