################################################
### 06-09-2022
##
#
# /mnt/extraids/OceanStor-0/linpei/hifi/data_6/lib1
L1R1=./raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
flowcell=000H3NWWM
i=1
for j in 1 2
do
surf=$flowcell:$i:$j
date
surfdedup $surf $L1R1 > L1R1Dedup_$i$j.fasta
date
finduniqread.pl L1R1Dedup_$i$j.fasta > L1R1Uniq_$i$j.fasta
done
###
L2R1=./raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
nohup bwa index -p L2R1 $L2R1 > bwaindexL2R1o 2>bwaindexL2R1e &
###
flowcell=000H3NWWM
i=1
# Mini-seq has one lane, on which there is top(j = 1) and bottom(j = 2) surfaces.
for j in 1 2
do
seq=L1R1Uniq_1$j
bwa mem L2R1 ../lib1/$seq.fasta -a -k $k -t 64 > $seq\_L2R1_ak$k.sam 2>$seq\_L2R1_ak$k.same
sam=$seq\_L2R1_ak$k.sam
grep -P "\t0\tMN00185:" $sam | cut -f 1,2,3 > $seq\_L2R1_ak$k\_mappedspot.L
grep -P "\t256\tMN00185:" $sam | cut -f 1,2,3 >> $seq\_L2R1_ak$k\_mappedspot.L
n1=`cat $seq\_L2R1_ak$k\_mappedspot.L | wc -l`
hifia_1 $seq\_L2R1_ak$k\_mappedspot.L $n1 > $seq\_L2R1_ak$k\_hifia_1.o
n2=`cat $seq\_L2R1_ak$k\_hifia_1.o | wc -l`
if [ -e hifia_2_$seq\_L2R1_ak$k\_mappedspot_out ]
then 
rm hifia_2_$seq\_L2R1_ak$k\_mappedspot_out
fi
for column in 1 2 3
do
for row in 101 102 103 104
do
tile=$flowcell:$i:$j$column$row
hifia_2 $seq\_L2R1_ak$k\_mappedspot.L $n1 $seq\_L2R1_ak$k\_hifia_1.o $n2 $tile >> hifia_2_$seq\_L2R1_ak$k\_mappedspot_out 2>>anye
done
done
done










#
##
##
################################################


run=10
L1R1=/mnt/extraids/OceanStor-0/linpei/hifi/data_$run/lib1/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
L2R1=/mnt/extraids/OceanStor-0/linpei/hifi/data_$run/lib2/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
L2R2=/mnt/extraids/OceanStor-0/linpei/hifi/data_$run/lib2/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz

# We had a 48 bps lab-made barcode
# ATVNNNNBATAVNNNNBTTAVNNNNBATTVNNNNBTATVNNNNVBNNN
# which should be present in HiFi R1 reads and used flow cell R1 reads.
# 
###
bwa index -p L1R1 $L1R1 > bwaindexo 2>bwaindexe
bwa mem L1R1 $L2R1 -a -t 64 > L2R1toL1R1_k0.sam 2>L2R1toL1_k0.e

grep -P "^MN00185:" L2R1toL1R1_k0.sam | cut -f 1,2,3 > L2R1toL1R1_k0.sam.column123.L

checkbarcode.pl L2R1.fastq > L2R1RAW_barcoded_check_result.txt
checkbarcode.pl L1R1.fastq > L1R1RAW_barcoded_check_result.txt

# R codes:
  sam = read.table("L2R1toL1R1_k0.sam.column123.L",sep="\t")
  spot2barcode = read.table("L1R1RAW_barcoded_check_result.txt",sep="\t");
  hifi2barcode = read.table("L2R1RAW_barcoded_check_result.txt",sep="\t")
  unique(sam[,2])
  for(i in unique(sam[,2])) {
    cat(i,sum(sam[,2] == i),"\n")
  }
  hifi2spot_1 = unique(sam[which(sam[,2] == 0 | sam[,2] == 256),1])
  hifi0spot_1 = unique(sam[which(sam[,2] == 4),1])
  sum(hifi0spot_1 %in% hifi2spot_1)
  sum(hifi2spot_1 %in% hifi2barcode[,1])/length(hifi2spot_1)
  sum(hifi0spot_1 %in% hifi2barcode[,1])/length(hifi0spot_1)
###
###
for k in 24 48
do
bwa mem L1R1 $L2R1 -a -k $k -t 64 > L2R1toL1_k$k.sam 2>L2R1toL1_k$k.e;
done

### count raw reads and deduplicated reads
##
# readedup: ~/bin
gunzip -c $L1R1 | grep "MN00185" | wc -l
gunzip -c $L2R1 | grep "MN00185" | wc -l

##
# examine the occurrence of 48 bps lab-made barcode in HiFi R1 reads and used flow cell R1 reads.
#
for i in 1 2
do
nohup perl readcheck_v2.pl barcode5 L$i\R1.fastq > L$i\R1barcodesurvey_v2_5.o 2>L$i\R1barcodesurvey_v2_5.e &
nohup perl readcheck_v2.pl barcode3 L$i\R1.fastq > L$i\R1barcodesurvey_v2_3.o 2>L$i\R1barcodesurvey_v2_3.e &
done

#
#
readedup $L1R1 > L1R1dedup.fasta 2>L1R1dedup.e
grep "MN00185" L1R1dedup.fasta | wc -l
readedup $L2R1 > L2R1dedup.fasta 2>L2R1dedup.e
grep "MN00185" L2R1dedup.fasta | wc -l
#
#
# map hifi reads (R2) to genes
#
# only used genes with gene name in gtf.
bwa mem $mwd/genome/release105/DNAMM39 $L2R2 -t 64 > bwaL2R2tomm39.sam 2>anye;
getgenefromgtf.pl $mwd/genome/release105/Mus_musculus.GRCm39.105.gtf ENSMUSG > genensmusg105.b 2>>anye
filetag=bwaL2R2tomm39
sam=$filetag.sam
bam=$filetag.bam
samtools view -S -b $sam --threads 16 > $bam 2>>anye
genicreadfile=$filetag\gene.L
bedtools intersect -a $bam -b genensmusg105.b -wb -bed > $genicreadfile
cut -f 1,2,3,4,16 $genicreadfile > $filetag\gene.fivecolumn.L

cut -f 4 $genicreadfile | sort | uniq > hifireads_biologically_resolved.L
cut -f 16 $genicreadfile | sort | uniq > hifireads_biologically_resolved_genes.L
grep "protein_coding" hifireads_biologically_resolved_genes.L | wc -l
grep "lncRNA" hifireads_biologically_resolved_genes.L | wc -l


#
###
##
# map HiFi R1 reads to used flow cell R1 reads using BWA MEM
for k in 24 48
do
bwa mem L1R1 $L2R1 -a -k $k -t 64 > L2R1toL1_k$k.sam
grep -P "\t0\tMN00185:" L2R1toL1_k$k.sam | cut -f 1,2,3 > L2R1toL1_k$k\_0256.cleansam 
grep -P "\t256\tMN00185:" L2R1toL1_k$k.sam | cut -f 1,2,3 >> L2R1toL1_k$k\_0256.cleansam 
cut -f 1 L2R1toL1_k$k\_0256.cleansam | sort | uniq -c > hifireads_spatially_resolved_n_spot_k$k.L
n1=`grep -P "\s1\sMN00185" hifireads_spatially_resolved_n_spot_k$k.L | wc -l`
n2=`grep -P "\s2\sMN00185" hifireads_spatially_resolved_n_spot_k$k.L | wc -l`
nn=`grep -P "\s\d+\sMN00185" hifireads_spatially_resolved_n_spot_k$k.L | wc -l`
let "n3=nn-n1-n2"
echo $k $nn $n1 $n2 $n3
done
# nn: number of hifi reads (R1) could be mapped to any spot on the used flow cell.
# n1: same as nn, but only hifi reads mapped to 1 spot.
# n2: 2 spot
# n3: >= 3 spot
##
###
###
##
#

L1R1=./lib1/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
if [ -e n_raw_spots_per_tile.L ]
then
rm n_raw_spots_per_tile.L
fi
for i in `cut -f 1 Tiles.L`
do
n=`gunzip -c $L1R1 | grep "MN00185:251:000H3VFVV:1:$i:" | wc -l`
echo $i $n >> n_raw_spots_per_tile.L
done

if [ -e n_hifi_per_tile.L ]
then
rm n_hifi_per_tile.L
fi

# i: the ID of each tile. 
# in this case,000H3VFVV is the ID of the used flow cell.
for i in `cut -f 1 Tiles.L`
do
  grep -P "\t0\tMN00185:251:000H3VFVV:1:$i:" L2R1toL1_k24.sam | cut -f 1 >  n_hifi_reads_spatially_resolved_per_Tile$i.L
grep -P "\t256\tMN00185:251:000H3VFVV:1:$i:" L2R1toL1_k24.sam | cut -f 1 >> n_hifi_reads_spatially_resolved_per_Tile$i.L
sort n_hifi_reads_spatially_resolved_per_Tile$i.L | uniq > hifi_reads_spatially_resolved_per_Tile$i.L
rm n_hifi_reads_spatially_resolved_per_Tile$i.L
n=`cat hifi_reads_spatially_resolved_per_Tile$i.L | wc -l`
echo $i $n >> n_hifi_per_tile.L
rm hifi_reads_spatially_resolved_per_Tile$i.L
done
#
##
###
if [ -e n_hifi_genic_per_tile.L ]
then
rm n_hifi_genic_per_tile.L
fi

# working dir:
# /mnt/extraids/OceanStor-0/linpei/hifi/data_10

# the follow counted the number of gene-mapped spatially resolved HiFi reads.
grep -P "\t0\tMN00185:251:000H3VFVV:1:" L2R1toL1_k24.sam | cut -f 1 > spatResolved_hifi.L
grep -P "\t256\tMN00185:251:000H3VFVV:1:" L2R1toL1_k24.sam | cut -f 1 >> spatResolved_hifi.L
sort spatResolved_hifi.L | uniq > spatResolved_hifi_unique.L
cut -f 4 bwaL2R2tomm39gene.fivecolumn.L | sort | uniq > geneMapped_hifi_uniq.L
cat spatResolved_hifi_unique.L geneMapped_hifi_uniq.L | sort | uniq -c | grep -P "2\sMN" | wc -l
rm spatResolved_hifi.L spatResolved_hifi_unique.L geneMapped_hifi_uniq.L

################################################
################################################ stringent: barcoded or uniquely spatially resolved
################################################
n2=`cat bwaL2R2tomm39_L2R2maskedgene.fivecolumn.L | wc -l`
k=48
for C in barcoded hifi1spot
do
n1=`cat L2R1toL1_k$k\_0256_$C.cleansam | wc -l`
if [ -e hifia2_v0513_k$k\_$C.o ]
then
rm hifia2_v0513_k$k\_$C.o
fi
for i in `cut -f 1 Tiles.L`; 
do 
tile=000H3VFVV:1:$i: 
./hifia2 L2R1toL1_k$k\_0256_$C.cleansam $n1 bwaL2R2tomm39_L2R2maskedgene.fivecolumn.L $n2 $tile >> hifia2_v0513_k$k\_$C.o
done
done
################################################



n1=`cat L2R1toL1_k24_0256.cleansam | wc -l`
n2=`cat bwaL2R2tomm39gene.fivecolumn.L | wc -l`
for i in `cut -f 1 Tiles.L`; 
do 
tile=000H3VFVV:1:$i: 
date; 
# the output in HiFiSlide_hifi_to_Tile$i\_by_gene.o has 4 columns:
# HiFi read, number of mapped spots on tile $i, spot read (from the used flow cell), gene
./hifia L2R1toL1_k24_0256.cleansam $n1 bwaL2R2tomm39gene.fivecolumn.L $n2 $tile > HiFiSlide_hifi_to_Tile$i\_by_gene.o 2>HiFiSlide_hifi_to_Tile$i\_by_gene.e
# reform the output to plot dots of genes on the tile.
# cut -f 3,4 HiFiSlide_hifi_to_Tile$i\_by_gene.cpp | perl -p -e "s/:/\t/g" | cut -f 1,2,3,5,6,7,8,13 > HiFiSlide_hifi_to_Tile$i\_by_gene.dot
#
nn=`grep -P "\t\d+\tMN00185" HiFiSlide_hifi_to_Tile$i\_by_gene.o | cut -f 1 | sort | uniq | wc -l`
echo $i $nn >> n_hifi_genic_per_tile.L
date; 
done
### source code: HiFianalysis.cpp in /mnt/extraids/OceanStor-0/linpei/hifi/data_10
## compile: g++ -o hifia HiFianalysis.cpp -lz
##
# working dir:
# /mnt/extraids/OceanStor-0/linpei/hifi/data_10
plot_Xnspot_Ynhifi_per_tile.R

# input for plot_Xnspot_Ynhifi_per_tile.R:
# (1) n_hifi_per_tile.L: number of hifi read pairs whose R1 could be mapped to any spot on each tile
# (2) n_raw_spots_per_tile.L: number of raw spots on each tile
# (3) n_hifi_genic_per_tile.L: same as (1) but only considered hifi read pairs whose R2 could be mapped to any genes.

#
# TO BE CONTINUED
#
# 
#
##############################################
###################################################################################################################
#
# gunzip -c $L1R1 | grep "MN00185" | cut -d " " -f 1 | cut -d ":" -f 5 | sort | uniq

nohup bwa index -p L1R1 $L1R1 > bwaindexo 2>bwaindexe &

k=24
bwa mem L1R1 $L2R1 -a -k $k -t 64 > L2R1toL1_k$k.sam 2>L2R1toL1_k$k.e;


grep -P "^MN00185:260:000H3W2LL:\d:\d+:\d+:\d+\t0\t" L2R1toL1_k$k.sam | cut -f 1,2,3 > L2R1toL1_k$k.hits
grep -P "^MN00185:260:000H3W2LL:\d:\d+:\d+:\d+\t256\t" L2R1toL1_k$k.sam | cut -f 1,2,3 >> L2R1toL1_k$k.hits

gethificoord.py L2R1toL1_k$k.hits > one_hit_hifi_reads_nspot2tiles.L 2>n_hits_hifi_reads.L
# sort -k 2nr one_hit_hifi_reads_nspot2tiles.L

################################################################################################################
date;
countreadundancy $L1R1 > Library1_readundancy_check.o 2>Library1_readundancy_check.e;
date
if [ -e unique_coord_per_tile.o ]; 
then 
echo "Found";
rm unique_coord_per_tile.o; 
fi

for j in `gunzip -c $L1R1 | grep "MN00185" | cut -d " " -f 1 | cut -d ":" -f 5 | sort | uniq`;
do
tile=FVV:1:$j:
n=`grep $tile Library1_readundancy_check.o | grep -P "\s1$" | wc -l` 
echo $tile $n >> unique_coord_per_tile.o
done
################################################################################################################
# 04-21
L2R1=./lib2/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
L2R2=./lib2/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz

# 1: how many barcode with "good" seq can be mapped by "good" hifi reads?
#
# we had a 48 bps long structure and it SHOULD be present
perl readcheck.pl L1R1.fastq > L1R1good48.L
perl readcheck.pl L2R1.fastq > L2R1good48.L

echo "#" > L1R1good48_n_of_barcode_per_tile.L
for i in `grep -P "\t1$" L1R1good48.L | cut -d ":" -f 5 | sort | uniq`;
do
tile=:1:$i:
n=`grep $tile L1R1good48.L | grep -P "\t1$" | wc -l`;
echo $i $tile $n >> L1R1good48_n_of_barcode_per_tile.L
done

k=35
bwa mem L1R1 $L2R1 -a -t 64 > L2R1toL1_k$k.sam 2>L2R1toL1_k$k.e;
grep -P "^MN00185:260:000H3W2LL:1:\d+:\d+:\d+\t0\t" L2R1toL1_k$k.sam | cut -f 1,2,3 > L2R1toL1_k$k.hits
grep -P "^MN00185:260:000H3W2LL:1:\d+:\d+:\d+\t256\t" L2R1toL1_k$k.sam | cut -f 1,2,3 >> L2R1toL1_k$k.hits
gethificoord.py L2R1toL1_k$k.hits L1R1good48.L L2R1good48.L > one_hit_hifi_reads_nspot2tiles_k$k.L 2>n_hits_hifi_reads_k$k.L

# 2: how many "good" hifi can be mapped to genes?



cut -f 4 $genicreadfile | sort | uniq > bwaL2R2tomm39genicread.L
cut -f 1 L2R1good48.L | sort | uniq > L2R1good48_uniq.L
cat L2R1good48_uniq.L > mug
cat bwaL2R2tomm39genicread.L >> mug
sort mug | uniq -c | grep -P "1\sMN00185" | wc -l
sort mug | uniq -c | grep -P "2\sMN00185" | wc -l
sort mug | uniq -c | grep -P "MN00185" | wc -l


#############################
gunzip -c $L2R1 > L2R1.fastq
gunzip -c $L1R1 > L1R1.fastq
head n_hits_hifi_reads_kdefault.L > templines
 date;getreadsbylist.py $L2R1 templines 2>anye;date
Wed Apr 20 12:30:27 PDT 2022
>MN00185:260:000H3W2LL:1:11101:15068:1055|1
GCATCAGTTGATAATGTGTTTANCCANTATTGNTNNTTNTAATGNATNNTGNNAN
>MN00185:260:000H3W2LL:1:11101:22712:1058|1
GCTTGTAAACATAGGTAAGTTAACGTNTATTGNAANTTNTAGTGGATCNGGNNAN
>MN00185:260:000H3W2LL:1:11101:5740:1061|1
GCCTAAGCACATAAATCATTTAGGCGNTATTGTATCTTNTGTCAAGTGNGGNGAN
>MN00185:260:000H3W2LL:1:11101:14859:1062|112
CCCCTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTNGAAAAN
>MN00185:260:000H3W2LL:1:11101:21044:1063|112
CCCCTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTNGACAAN
>MN00185:260:000H3W2LL:1:11101:8269:1064|112
CCCCTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTNGAAAAN
>MN00185:260:000H3W2LL:1:11101:7972:1065|1
CCCCTTCAGCAGACATGCCGAGCCGTATCTCGTATGCCGTCTTCTGCTNGAAAAN
>MN00185:260:000H3W2LL:1:11101:19559:1065|1
AGAGTGACCAGACATGCCGAGCCGTATCTCGTATGCCGTCTTCTGCTTNAAAAAN
>MN00185:260:000H3W2LL:1:11101:17366:1066|1
CCCCTTCACCAGACATGCCGAGCCGTATCTCGTATGCCGTCTTCTGCTNGAAAAN
>MN00185:260:000H3W2LL:1:11101:10811:1066|1
GCAGTGACCAGACATGCCGAGCCGTATCTCGTATGCCGTCTTCTGCTTNAAAATN
Wed Apr 20 12:30:54 PDT 2022

grep -P "^MN00185:260:000H3W2LL:1:11101:14859:1062" L2R1toL1_k.sam > checkfile1
cut -f 3 checkfile1 > temolines2
date;getreadsbylist.py L1R1.fastq temolines2 > L1R1templines2.fasta 2>anye;date


    GCATGCCATGATACGATTGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAA
    GCATGGAATGATAGGGGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAA
        GCATCATGATATACGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAA
            GCAGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAGG
     GCATCATTGTATACTATTGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAA

HiFi                CCCCTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTNGAAAAN
