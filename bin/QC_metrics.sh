###### TO BE APPENDED TO THE END OF HiFi_processing_pipeline.sh

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
a=$(wc -l $L1_DIR/L1R1_dedup.fasta | cut -d " " -f 1)
m2=$(($a / 2))

a=$(echo "scale=4 ; $m2 / $m1 * 100" | bc | awk '{printf("%.2f",$1)}')
m3=$a"%"

##### Number of input HiFi read pairs
m4=$(grep -w "M\\:\\:mem_process_seqs" $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.log |
cut -d " " -f3 | xargs | tr ' ' + | bc)

##### Number of HiFi read pairs aligned to spatial barcodes (number of unique HiFi read IDs in column 1 of L2R1_L1R1_dedup.hifislida.o)
m5=$(awk -F '\t' '{print $1}' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.o | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m5 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m6=$a"%"


##### Number of HiFi read pairs aligned to unique spatial barcodes (number of HiFi read IDs in column 1 of L2R1_L1R1_dedup.hifislida.o with column 3 and 4 equal to 1)
m7=$(cut -f3 $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o | xargs | tr ' ' + | bc) # same as below but faster
# awk -F "\t" '$3 == 1 && $4 == 1 { print $1 }' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup_1_1.hifislida.o | wc -l

a=$(echo "scale=4 ; $m7 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m8=$a"%"

##### Number of HiFi read pairs aligned to unique spatial barcodes and under ROI
# awk 'BEGIN{FS=OFS="\t"} {gsub(/T/, "", $2)} 1' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup_1_1.hifislida2.o | head -10 # to remove the T (not used anymore)
m9=$(awk -F "\t" 'NR==FNR{a[$1]; next} FNR==0 || $2 in a' $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o | cut -f3 | xargs | tr ' ' + | bc)

a=$(echo "scale=4 ; $m9 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m10=$a"%"


##### Number of HiFi read pairs under ROI (all spatially resolved, i.e. not only aligned to unique spatial barcodes)
m11=$(awk -F '\t' 'NR>1 {print $1}' $L2_DIR/L2R1_mapping/L2R1_L1R1.hifislida3.o | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m11 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m12=$a"%"


##### Number of HiFi-Slide L2R2 passing pear filtering
a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2_pear.unassembled.reverse.fastq | cut -d " " -f 1)
m13=$(($a / 4))

a=$(echo "scale=4 ; $m13 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m14=$a"%"

##### Number of HiFi-Slide L2R2 passing length filtering (performed automatically by fastp)
# a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq | cut -d " " -f 1)
# echo $(($a / 4))
# To make it faster, we can just count the number of input reads in the log of the STAR aligner
m15=$(grep -w "Number of input reads" $L2_DIR/L2R2_mapping/genome/L2R2_genome.Log.final.out |cut -d "|" -f2 | sed 's/\t//g')

a=$(echo "scale=4 ; $m15 / $m4 * 100" | bc | awk '{printf("%.2f",$1)}')
m16=$a"%"


########## GENOME

##### Number of HiFi-Slide L2R2 uniquely mapped to genome
# grep -w "Uniquely mapped reads number" $L2_DIR/L2R2_mapping/genome/L2R2_genome.Log.final.out | cut -d "|" -f2 | sed 's/\t//g'

##### Number of HiFi-Slide L2R2 uniquely mapped to genome and to annotated genes
m17=$(cut -f4 $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m17 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m18=$a"%"

##### Number of HiFi-Slide read pairs genome mapped and spatially resolved
m19=$(awk 'FNR==NR {a[$1]; next} FNR> 0 && $4 in a' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.o $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed | cut -f4 | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m19 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m20=$a"%"


##### Number of HiFi reads spatially resolved, under ROI and aligned to genome
m21=$(awk -F '\t' 'NR>1 {print $1}' $L2_DIR/L2R1_L2R2_integrate/HiFi_L2R2_genome_spatial.txt | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m19 / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m22=$a"%"

# Number of genes 
m23=$(awk -F '\t' 'NR>1 {print $9}' $L2_DIR/L2R1_L2R2_integrate/HiFi_L2R2_genome_spatial.txt | sort --parallel=$N_THREADS | uniq | wc -l)

##### Average number of HiFi-Slide read pairs genome mapped and spatially resolved per tile
n_tiles=$(wc -l $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o | cut -d " " -f 1)
m24=$(echo "scale=4 ; $m19 / $n_tiles" | bc | awk '{printf("%.2f",$1)}')

##### Average number of HiFi-Slide read pairs genome mapped and spatially resolved per tile under ROI
n_tiles_under_ROI=$(wc -l $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt | cut -d " " -f 1)
m25=$(echo "scale=4 ; $m21 / $n_tiles_under_ROI" | bc | awk '{printf("%.2f",$1)}')


#### Labels
M1="Total number of barcodes"
M2="Number of deduplicated barcodes"
M3="% of deduplicated barcodes"
M4="Number of HiFi-Slide read pairs"
M5="Number of HiFi-Slide read pairs spatially resolved (aligned to spatial barcodes)"
M6="% of HiFi-Slide read pairs spatially resolved"
M7="Number of HiFi-Slide read pairs univocally spatially resolved (aligned to unique spatial barcodes)"
M8="% of HiFi-Slide read pairs univocally spatially resolved"
M9="Number of HiFi-Slide read pairs univocally spatially resolved and under ROI"
M10="% of HiFi-Slide read pairs univocally spatially resolved and under ROI"
M11="Number of HiFi-Slide read pairs spatially resolved (aligned to spatial barcodes) and under ROI"
M12="% of HiFi-Slide read pairs spatially resolved (aligned to spatial barcodes) and under ROI"
M13="Number of HiFi-Slide read pairs passing PEAR filtering"
M14="% of HiFi-Slide read pairs passing PEAR filtering"
M15="Number of HiFi-Slide read pairs passing PEAR and FASTP (length) filtering"
M16="% of HiFi-Slide read pairs passing PEAR and FASTP (length) filtering"
M17="Number of HiFi-Slide read pairs genome mapped (uniquely aligned to annotated genes)"
M18="% of HiFi-Slide read pairs genome mapped"
M19="Number of HiFi-Slide read pairs genome mapped and spatially resolved"
M20="% of HiFi-Slide read pairs genome mapped and spatially resolved"
M21="Number of HiFi-Slide read pairs genome mapped and under ROI"
M22="% of HiFi-Slide read pairs genome mapped and under ROI"
M23="Number of annotated genes aligned to HiFi-Slide read pairs"
M24="Average number of HiFi-Slide read pairs genome mapped and spatially resolved per tile"
M25="Average number of HiFi-Slide read pairs genome mapped and spatially resolved per tile under ROI"


rm $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
touch $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
for k in $(seq 1 25); do
Mk=M${k}
mk=m${k}
echo -e ${!Mk}'\t'${!mk}>> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
done

# printf '%s\n' "$M1" "$M2" "$M3" "$M4" "$M5" "$M6" "$M7" "$M8" "$M9" "$M10" "$M11" "$M12" "$M13" "$M14" "$M15" "$M16" "$M17" "$M18" "$M19" "$M20" "$M21" | paste -sd '\n' > $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log




########## TRANSCRIPTOME
for my_transcript in tRNA piRNA miRNA circRNA; do

size=$(stat -c %s $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt")

if [ $size != 0 ]; then

##### Number of uniquely mapped reads to the transcriptome
m26a=$(grep -w "aligned exactly 1 time" $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_mapped.log" | cut -d " " -f5)

a=$(echo "scale=4 ; $m26a / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m26b=$a"%"

##### Number of HiFi-Slide read pairs transcriptome mapped and spatially resolved
m26c=$(awk 'FNR==NR {a[$1]; next} FNR> 0 && $1 in a' $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida.o $L2_DIR/L2R2_mapping/transcriptome/$my_transcript/L2R2_$my_transcript"_uniquely_mapped.txt" | cut -f1 | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m26c / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m26d=$a"%"

##### Number of HiFi reads spatially resolved, under ROI and aligned to transcriptome
m26e=$(awk -F '\t' 'NR>1 {print $1}' $L2_DIR/L2R1_L2R2_integrate/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" | sort --parallel=$N_THREADS | uniq | wc -l)

a=$(echo "scale=4 ; $m26c / $m15 * 100" | bc | awk '{printf("%.2f",$1)}')
m26f=$a"%"

# Number of transcripts 
m26g=$(awk -F '\t' 'NR>1 {print $6}' $L2_DIR/L2R1_L2R2_integrate/$my_transcript/HiFi_L2R2_$my_transcript"_spatial.txt" | sort --parallel=$N_THREADS | uniq | wc -l)

##### Average number of HiFi-Slide read pairs transcriptome mapped and spatially resolved per tile
n_tiles=$(wc -l $L2_DIR/L2R1_mapping/L2R1_L1R1_dedup.hifislida2.o | cut -d " " -f 1)
m26h=$(echo "scale=4 ; $m26c / $n_tiles" | bc | awk '{printf("%.2f",$1)}')

##### Average number of HiFi-Slide read pairs transcriptome mapped and spatially resolved per tile under ROI
n_tiles_under_ROI=$(wc -l $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt | cut -d " " -f 1)
m26i=$(echo "scale=4 ; $m26e / $n_tiles_under_ROI" | bc | awk '{printf("%.2f",$1)}')

elif [ $size == 0 ]; then
m26a=0
m26b="0%"
m26c=0
m26d="0%"
m26e=0
m26f="0%"
m26g=0
m26h=0
m26i=0

fi

M26a=$my_transcript" - Number of HiFi-Slide read pairs transcriptome mapped (uniquely aligned to transcripts)"
M26b=$my_transcript" - % of HiFi-Slide read pairs transcriptome mapped"
M26c=$my_transcript" - Number of HiFi-Slide read pairs transcriptome mapped and spatially resolved"
M26d=$my_transcript" - % of HiFi-Slide read pairs transcriptome mapped and spatially resolved"
M26e=$my_transcript" - Number of HiFi-Slide read pairs transcriptome mapped and under ROI"
M26f=$my_transcript" - % of HiFi-Slide read pairs transcriptome mapped and under ROI"
M26g=$my_transcript" - Number of transcripts aligned to HiFi-Slide read pairs"
M26h="Average number of HiFi-Slide read pairs transcriptome mapped and spatially resolved per tile"
M26i="Average number of HiFi-Slide read pairs transcriptome mapped and spatially resolved per tile under ROI"

for k in a b c d e; do
Mk=M26${k}
mk=m26${k}
echo -e ${!Mk}'\t'${!mk}>> $OUT_DIR/$SAMPLE_NAME/$SAMPLE_NAME.log
done

done
















