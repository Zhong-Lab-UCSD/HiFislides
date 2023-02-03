#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    ---------------------------------------------------------------------------
    Usage: $PROGNAME [-R <roi_tiles>] [-r <roi_label>] [-N <sample_name>]
                     [-o <output_dir>]
    
    This is a script to extract the tiles under ROI from the final outpur file of the HiFi Pipeline.
    -R : Path to the txt file with the tiles under ROI (each tile on a new line).
    -r : ROI name to label the output file.
    -N : Sample name, used to label the final output file.
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

while getopts :R:r:N:o:h opt; do
    case $opt in
        R) ROI_TILES=${OPTARG};;
        r) ROI_label=${OPTARG};;
        N) SAMPLE_NAME=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        h) usage;;
    esac
done

# check parameters
[ -z "$ROI_TILES" ] &&  echo "Error!! Please provide path to the txt file with the tiles under ROI (each tile on a new line) with -R" && parameter_error

[ -z "$ROI_label" ] && echo "Error!! Please provide the ROI label to be used with -r" && parameter_error

[  -z "$SAMPLE_NAME" ] && echo "Error!! Please provide the sample name with -N" && parameter_error

[ ! -d "$OUT_DIR" ] && echo "Error!! Output directory not exist: "$OUT_DIR && parameter_error


##############################
L2R1_L2R2_INTEGRATE_DIR=$OUT_DIR/L2R1_L2R2_integrate

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1"_"$2"_"$3, $1, $2, $3, $5, $6, $7}' | awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]+=$7} END { for (i in array) {print i"\t" array[i]}}' > $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.$ROI_label.txt




















