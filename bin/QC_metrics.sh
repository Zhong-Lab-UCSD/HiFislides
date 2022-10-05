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

##### Number of HiFi read pairs  aligned to unique spatial barcodes and under ROI
# awk 'BEGIN{FS=OFS="\t"} {gsub(/T/, "", $2)} 1' $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup_1_1.hifislida2.o | head -10 # to remove the T (not used anymore)
awk -F "\t" 'NR==FNR{a[$1]; next} FNR==0 || $2 in a' $L2_DIR/L2R1_mapping/ROI_tile_IDs.txt $L2_DIR/L2R1_mapping/L2R1__L1R1_dedup.hifislida2.o | cut -f3 | xargs | tr ' ' + | bc


##### Number of HiFi read pairs under ROI (all spatially resolved, i.e. not only aligned to unique spatial barcodes)
awk -F '\t' 'NR>1 {print $1}' $L2_DIR/L2R1_mapping/L2R1__L1R1.hifislida3.o | sort --parallel=32 | uniq | wc -l


##### Number of HiFi-Slide L2R2 passing pear filtering
a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2_pear.unassembled.reverse.fastq | cut -d " " -f 1)
echo $(($a / 4))

##### Number of HiFi-Slide L2R2 passing length filtering (performed automatically by fastp)
# a=$(wc -l $L2_DIR/L2R2_preprocessing/L2R2.trim_front_60.fastq | cut -d " " -f 1)
# echo $(($a / 4))
# To make it faster, we can just count the number of input reads in the log of the STAR aligner
grep -w "Number of input reads" $L2_DIR/L2R2_mapping/genome/L2R2_genome.Log.final.out |cut -d "|" -f2 | sed 's/\t//g'

########## GENOME

##### Number of HiFi-Slide L2R2 uniquely mapped to genome
grep -w "Uniquely mapped reads number" $L2_DIR/L2R2_mapping/genome/L2R2_genome.Log.final.out | cut -d "|" -f2 | sed 's/\t//g'

##### Number of HiFi-Slide L2R2 uniquely mapped to genome and to annotated genes
wc -l $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed | cut -d " " -f 1

##### Number of HiFi reads spatially resolved, under ROI and aligned to genome
awk -F '\t' 'NR>1 {print $1}' $L2_DIR/L2R1_L2R2_integrate/HiFi_L2R2_genome_spatial.txt | sort --parallel=32 | uniq | wc -l

# Number of genes 
awk -F '\t' 'NR>1 {print $9}' $L2_DIR/L2R1_L2R2_integrate/HiFi_L2R2_genome_spatial.txt | sort --parallel=32 | uniq | wc -l


########## TRANSCRIPTOME

for t in tRNA piRNA mirbase circbase; do

##### Number of uniquely mapped reads to the transcriptome
grep -w "aligned exactly 1 time" $L2_DIR/L2R2_mapping/transcriptome/$t/L2R2_$t"_mapped.log" | cut -d " " -f5

##### Number of HiFi reads spatially resolved, under ROI and aligned to transcriptome
awk -F '\t' 'NR>1 {print $1}' $L2_DIR/L2R1_L2R2_integrate/$t/HiFi_L2R2_$t"_spatial.txt" | sort --parallel=32 | uniq | wc -l

# Number of transcripts 
awk -F '\t' 'NR>1 {print $6}' $L2_DIR/L2R1_L2R2_integrate/$t/HiFi_L2R2_$t"_spatial.txt" | sort --parallel=32 | uniq | wc -l

done













