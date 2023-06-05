#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    ---------------------------------------------------------------------------
    Usage: $PROGNAME [-R <roi_tiles>] [-r <roi_label>] [-N <sample_name>]
                     [-i <integrate_dir>]
    
    This is a script to extract the tiles under ROI from the final outpur file of the HiFi Pipeline.
    -R : Path to the txt file with the tiles under ROI (each tile on a new line).
    -r : ROI name to label the output file.
    -N : Sample name, used to label the final output file.
    -i : Input directory where integrated data are.
    -h : Show usage help
    ---------------------------------------------------------------------------
EOF
    exit 1
}

parameter_error() {
    echo "!!!!!!!!!!!!!!!! Pipeline exit with parameter ERROR at $(date '+%Y-%m-%d %H:%M:%S %Z') !!!!!!!!!!!!!!!!"
    usage
}

while getopts :R:r:N:i:h opt; do
    case $opt in
        R) ROI_TILES=${OPTARG};;
        r) ROI_label=${OPTARG};;
        N) SAMPLE_NAME=${OPTARG};;
        i) INTEGRATE_DIR=${OPTARG};;
        h) usage;;
    esac
done

# check parameters
[ -z "$ROI_TILES" ] &&  echo "Error!! Please provide path to the txt file with the tiles under ROI (each tile on a new line) with -R" && parameter_error

[ -z "$ROI_label" ] && echo "Error!! Please provide the ROI label to be used with -r" && parameter_error

[  -z "$SAMPLE_NAME" ] && echo "Error!! Please provide the sample name with -N" && parameter_error

[ ! -d "$INTEGRATE_DIR" ] && echo "Error!! Input directory not exist: "$INTEGRATE_DIR && parameter_erro


##############################
echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] Start ROI extraction..."

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $ROI_TILES $INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1"_"$2"_"$3, $1, $2, $3, $5, $6, $7}' | awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]+=$7} END { for (i in array) {print i"\t" array[i]}}' > $INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.$ROI_label.txt

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d-%y %H:%M:%S')] ROI extraction finished."




















