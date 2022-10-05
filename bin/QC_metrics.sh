##### Total number of barcodes (L1R1)
touch $L1_DIR/L1R1_stats.txt

for x in $L1_FASTQ_DIR/$L1_FASTQ_BASENAME; do 
a=$(unpigz -p 32 -c $x | wc -l)
out=$(($a / 4))
echo -e $x"\t"$out >> $L1_DIR/L1R1_stats.txt
done

awk '{ sum += $2 } END { print sum }' $L1_DIR/L1R1_stats.txt

##### Number of deduplicated barcodes
a=$(wc -l $L1_DIR/L1R1_dedup.fasta | cut -d " " -f 1)
echo $(($a / 2))

##### Number of input HiFi read pairs
grep -w "M\\:\\:mem_process_seqs" $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup.log |
cut -d " " -f3 | xargs | tr ' ' + | bc

##### Number of HiFi read pairs aligned to spatial barcodes (number of unique HiFi read IDs in column 1 of L2R1__L1R1_dedup.hifislida.o)
awk -F '\t' '{print $1}' $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup.hifislida.o | sort --parallel=32 | uniq | wc -l

##### Number of HiFi read pairs aligned to unique spatial barcodes (number of HiFi read IDs in column 1 of L2R1__L1R1_dedup.hifislida.o with column 3 and 4 equal to 1)
cut -f3 $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_1_1.hifislida2.o | xargs | tr ' ' + | bc # same as below but faster
# awk -F "\t" '$3 == 1 && $4 == 1 { print $1 }' $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_1_1.hifislida.o | wc -l

##### Number of HiFi read pairs under ROI
# awk 'BEGIN{FS=OFS="\t"} {gsub(/T/, "", $2)} 1' $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_1_1.hifislida2.o | head -10 # to remove the T (not used anymore)
awk -F "\t" 'NR==FNR{a[$1]; next} FNR==0 || $2 in a' $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup.hifislida2.o | cut -f3 | xargs | tr ' ' + | bc

