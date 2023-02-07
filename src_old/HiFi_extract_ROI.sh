########## Select spots in ROI after visual inspection of the image
ROI_TILES=$L2_DIR/ROI_1_tiles.txt
ROI_label="ROI_1"

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1"_"$2"_"$3, $1, $2, $3, $5, $6, $7}' | awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]+=$7} END { for (i in array) {print i"\t" array[i]}}' > $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.$ROI_label.txt

ROI_TILES=$L2_DIR/ROI_2_tiles.txt
ROI_label="ROI_2"

awk -F"\t" 'NR==FNR{a[$1]; next} FNR==0 || $1 in a' $ROI_TILES $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.txt | awk -v OFS='\t' '{print $1"_"$2"_"$3, $1, $2, $3, $5, $6, $7}' | awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]+=$7} END { for (i in array) {print i"\t" array[i]}}' > $L2R1_L2R2_INTEGRATE_DIR/$SAMPLE_NAME.L2R2_genome_spatial.final.$ROI_label.txt



