#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    ---------------------------------------------------------------------------
    Usage: $PROGNAME [-N <sample_name>]
                     [-o <output_dir>] 
    
    This is a script to calculate HiFi stats.
    -N : Sample name.
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

while getopts :N:o:h opt; do
    case $opt in
        N) SAMPLE_NAME=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        h) usage;;
    esac
done

# check parameters
[  -z "$SAMPLE_NAME" ] && echo "Error!! Please provide the sample name with -N" && parameter_error

[ ! -d "$OUT_DIR" ] && echo "Error!! Output directory not exist: "$OUT_DIR && parameter_error


################## INITIALIZE VARIABLES

L2_DIR=$OUT_DIR/$SAMPLE_NAME # HiFi library
mkdir -p $L2_DIR

L2R2_GENOME_DIR=$L2_DIR/L2R2_mapping
mkdir -p $L2R2_GENOME_DIR

L2R1_MAPPING_DIR=$L2_DIR/L2R1_mapping
mkdir -p $L2R1_MAPPING_DIR

L2R1_L2R2_INTEGRATE_DIR=$L2_DIR/L2R1_L2R2_integrate
mkdir -p $L2R1_L2R2_INTEGRATE_DIR


################## STATS
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start stats calculation..."

##### Number of input HiFi read pairs
m1=$(grep -w "M\\:\\:mem_process_seqs" $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.log | cut -d " " -f3 | xargs | tr ' ' + | bc)

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
# m8=$(awk '!seen[$1]++' $L2R1_MAPPING_DIR/L2R1_L1R1_dedup.sam | wc -l)

# a=$(echo "scale=4 ; $m8 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
# m9=$a"%"

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
# M8="Number of read pairs spatially resolved"
# M9="Percentage of read pairs spatially resolved"
M10="Number of read pairs genome mapped and spatially resolved"
M11="Percentage of read pairs genome mapped and spatially resolved"
M12="Average spots per tile"
M13="Average spots per 10 um^2"
M14="Average genes per tile"
M15="Average genes per 10 um^2"

touch $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".stats.txt"
for k in $(1 2 3 4 5 6 7 10 11 12 13 14 15); do
# for k in $(seq 1 15); do
Mk=M${k}
mk=m${k}
echo -e ${!Mk}'\t'${!mk}>> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME".stats.txt"
done


echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] stats calculation finished." 
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] stats calculation finished." >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log


echo "Data processing started: "$START_DATE
echo "Data processing ended: "$(date)

echo "------------------------------" >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo "Data processing started: "$START_DATE >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
echo "Data processing ended: "$(date) >> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log



