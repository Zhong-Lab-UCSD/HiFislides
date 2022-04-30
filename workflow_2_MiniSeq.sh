
L1R1=./lib1/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
L2R1=./lib2/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
L2R2=./lib2/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz

#
gunzip -c $L1R1 | grep "MN00185" | wc -l
readedup $L1R1 > L1R1dedup.fasta 2>L1R1dedup.e
grep "MN00185" L1R1dedup.fasta | wc -l

gunzip -c $L2R1 | grep "MN00185" | wc -l
readedup $L2R1 > L2R1dedup.fasta 2>L2R1dedup.e
grep "MN00185" L2R1dedup.fasta | wc -l

#


bwa mem $mwd/genome/release105/DNAMM39 $L2R2 -t 64 > bwaL2R2tomm39.sam 2>anye;
getgenefromgtf.pl $mwd/genome/release105/Mus_musculus.GRCm39.105.gtf ENSMUSG > genensmusg105.b 2>>anye
filetag=bwaL2R2tomm39
sam=$filetag.sam
bam=$filetag.bam
samtools view -S -b $sam --threads 16 > $bam 2>>anye
genicreadfile=$filetag\gene.L
bedtools intersect -a $bam -b genensmusg105.b -wb -bed > $genicreadfile

cut -f 4 bwaL2R2tomm39gene.L | sort | uniq > hifireads_biologically_resolved.L
cut -f 16 bwaL2R2tomm39gene.L | sort | uniq > hifireads_biologically_resolved_genes.L
grep "protein_coding" hifireads_biologically_resolved_genes.L | wc -l
grep "lncRNA" hifireads_biologically_resolved_genes.L | wc -l

##
#
perl readcheck.pl L1R1.fastq > L1R1barcodesurvey.o
perl readcheck.pl L2R1.fastq > L2R1barcodesurvey.o

###
##
#
for k in 24 48
do 
grep -P "\t0\tMN00185:" L2R1toL1_k$k.sam | cut -f 1,2,3 > L2R1toL1_k$k\_0256.cleansam 
grep -P "\t256\tMN00185:" L2R1toL1_k$k.sam | cut -f 1,2,3 >> L2R1toL1_k$k\_0256.cleansam 
cut -f 1 L2R1toL1_k$k\_0256.cleansam | sort | uniq -c > hifireads_spatially_resolved_n_spot_k$k.L
n1=`grep -P "\s1\sMN00185" hifireads_spatially_resolved_n_spot_k$k.L | wc -l`
n2=`grep -P "\s2\sMN00185" hifireads_spatially_resolved_n_spot_k$k.L | wc -l`
nn=`grep -P "\s\d+\sMN00185" hifireads_spatially_resolved_n_spot_k$k.L | wc -l`
let "n3=nn-n1-n2"
echo $k $n1 $n2 $n3 $nn
done
#
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

