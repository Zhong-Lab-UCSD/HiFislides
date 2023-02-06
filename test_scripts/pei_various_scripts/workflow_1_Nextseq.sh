#
#################################################################
# /mnt/extraids/OceanStor-0/linpei/hifi/data_12/lib2
# run this as background process
#
flowcell=AAAHT3CHV
#
for i in 1 2
do
for j in 1 2
do
surf=$flowcell:$i:$j
date
surfdedup $surf *_L00$i\_R1_001.fastq.gz > L1R1Dedup_$i$j.fasta
date
finduniqread.pl L1R1Dedup_$i$j.fasta > L1R1Uniq_$i$j.fasta
done
done

for tile in 1314 1414 1514 1614 
do
finduniqread.pl L1R1Dedup_$i$j.fasta $tile > L1R1Uniq_$i$j\T$tile.fasta
echo $tile
done

#
# hifia_2_L1R1Uniq_11_L2R1_ak70_mappedspot_out
#################################################################

L2R1=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
L2R2=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz


date
# Tue May 31 11:56:44 PDT 2022
bwa index -p L2R1 $L2R1 > bwaindexL2R1o 2>bwaindexL2R1e
# Tue May 31 12:02:05 PDT 2022

k=40
seq=L1R1Uniq_11
sam=$seq\_L2R1_ak$k.sam
bwa mem L2R1 ../lib1/$seq.fasta -a -k $k -t 64 > $sam 2>$seq\_L2R1_ak$k.same

grep -P "\t0\tMN00185:" $sam | cut -f 1,2,3 > $seq\_L2R1_ak$k\_mappedspot.L
grep -P "\t256\tMN00185:" $sam | cut -f 1,2,3 >> $seq\_L2R1_ak$k\_mappedspot.L

###########
prefix=L2R2_000_

hg38=$mwd/imc/HG38
STAR --runThreadN 32 --genomeDir $hg38 --readFilesIn $L2R2 --quantMode GeneCounts --readFilesCommand zcat \ 
--outFileNamePrefix $prefix \ 
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 > starlogo 2>starlogoe
 
date
grep "SN:" $prefix\Aligned.out.sam > $prefix\Aligned.NH1.sam
grep ":STAR" $prefix\Aligned.out.sam >> $prefix\Aligned.NH1.sam 
grep -P "NH:i:1\t" $prefix\Aligned.out.sam >> $prefix\Aligned.NH1.sam 
date

filetag=$prefix\Aligned.NH1
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
gtf=$mwd/genome/release104/Homo_sapiens.GRCh38.104.gtf
getgenefromgtf.pl $gtf ENSG > genensg104.b 2>>anye
cat genensg104.b | perl -p -e "s/:/\t/g" | cut -f 1,2,3,4 > genensg104clean.b
bedtools intersect -a $bam -b genensg104clean.b -wb -bed > $genicreadfile
cut -f 1,2,3,4,16 $genicreadfile > $filetag\gene5columns.L

hifi2gene=L2R2_000_Aligned.NH1gene5columns.L
flowcell=AAAL33WM5
seq=L1R1Uniq_11
k=40

sam=$seq\_L2R1_ak$k.sam
hifia_asort.pl $sam $flowcell > $seq\_L2R1_ak$k\_mappedspot_1n.L;


hifia_1n_marker_per_spot.pl $seq\_L2R1_ak$k\_mappedspot_1n.L $hifi2gene $flowcell cell_type_marker_$i.txt $ensgname > Output_spot_to_gene_for_$i.txt
hifia_1n_marker_per_spot.pl $seq\_L2R1_ak$k\_mappedspot_1n.L $hifi2gene $flowcell hifi2gene.G $ensgname > Output_spot_to_gene.A

###########



nohup hifia_nmappedbc.pl $seq\_L2R1_ak$k\_mappedspot.L $flowcell > bc_per_tile.o 2>bc_per_tile.e &
nohup hifia_nmappedhifi.pl $seq\_L2R1_ak$k\_mappedspot.L $flowcell > hifi_per_tile.o 2>hifi_per_tile.e &

nohup hifia_ngenomichifir2_per_tile.pl $seq\_L2R1_ak$k\_mappedspot.L L2R2_000_Aligned.NH1.sam $flowcell > hifia_ngenomichifir2_per_tile.o &
nohup hifia_ngenichifir2_per_tile.pl $seq\_L2R1_ak$k\_mappedspot.L L2R2_000_Aligned.NH1gene5columns.L $flowcell > hifia_ngenichifir2_per_tile.o &













for tile in 1314 1414 1514 1614
do
seq=L1R1Uniq_11T$tile
bwa mem L2R1 $seq.fasta -a -k $k -t 64 > $seq\_L2R1_ak$k.sam 2>$seq\_L2R1_ak$k.same
hifi2bestalignedspot.pl $seq\_L2R1_ak$k.sam L2R2_010_Aligned.NH1.2column $tile > hifi2coord_$tile.o1 2>hifi2coord_$tile.o2
n=`cat hifi2coord_$tile.o1 | wc -l`
calcdist hifi2coord_$tile.o1 $n > hifi2dist_$tile.o
done











# ~/bin
# HiFianalysis_nmappedbarcode_per_tile.cpp
# hifianmbc
# HiFianalysis_ngenicbarcode_per_tile.cpp
# hifiangbc
# hifiangr2.pl: to reproduce data_8 results.
# input: 
# n1=`cat $seq\_L2R1_ak$k\_mappedspot.L | wc -l`
# $seq\_L2R1_ak$k\_mappedspot.L $n1
# n2=`cat L2R2_010_Aligned.NH1gene5columns.L | wc -l`
# L2R2_010_Aligned.NH1gene5columns.L

# cmd:
# out=hifi2gene2tile_06202022.o
# hifiangr2.pl L1R1Uniq_11_L2R1_ak70_mappedspot.L L2R2_015_Aligned.NH1gene5columns.L $out
# R CMD BATCH --no-save --no-restore "--args hifi2gene2tile_06202022.o 8 1000 64" testGeneToTile_1.R
# /mnt/extraids/OceanStor-0/linpei/hifi/data_8/lib2/plot_pval_per_tile.R


k=70
for seq in L1R1Uniq_11 L1R1Uniq_12 L1R1Uniq_21 L1R1Uniq_22
do
bwa mem L2R1 ../lib1/$seq.fasta -a -k $k -t 64 > $seq\_L2R1_ak$k.sam 2>$seq\_L2R1_ak$k.same
tile=1314
seq=L1R1Uniq_11T$tile
bwa mem L2R1 $seq.fasta -a -k $k -t 64 > $seq\_L2R1_ak$k.sam 2>$seq\_L2R1_ak$k.same
hifi2bestalignedspot.pl $seq\_L2R1_ak$k.sam L2R2_010_Aligned.NH1.2column $tile > $seq\_L2R1_ak$k_spot_per_hifi.o1 2>$seq\_L2R1_ak$k_spot_per_hifi.o2

sam=$seq\_L2R1_ak$k.sam
grep -P "\t0\tMN00185:" $sam | cut -f 1,2,3 > $seq\_L2R1_ak$k\_mappedspot.L
grep -P "\t256\tMN00185:" $sam | cut -f 1,2,3 >> $seq\_L2R1_ak$k\_mappedspot.L
done

# grep -v "ENSG00000230876" L2R2_010_Aligned.NH1gene5columns.L > G
for i in 1 2 3
do
for j in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do
nohup hifiangbc L1R1Uniq_11_L2R1_ak70_mappedspot.L $n1 AAAHT3CHV:1:1$i$j G $n2 > genicspotTile1$i$j.o 2>>anye &
done
done

if [ -e n_genic_spot_per_tile.o ]
then
rm n_genic_spot_per_tile.o
fi

for i in 1 2 3 4 5 6
do
for j in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do
tile=1$i$j
n=`cut -f 2 genicspotTile1$i$j.o | sort | uniq | wc -l`
echo $tile $n >> n_genic_spot_per_tile.o
done
done

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

################### Read 1
k=70
for seq in L1R1Uniq_11 L1R1Uniq_12 L1R1Uniq_21 L1R1Uniq_22
do
bwa mem L2R1 ../lib1/$seq.fasta -a -k $k -t 64 > $seq\_L2R1_ak$k.sam 2>$seq\_L2R1_ak$k.same
sam=$seq\_L2R1_ak$k.sam
grep -P "\t0\tMN00185:" $sam | cut -f 1,2,3 > $seq\_L2R1_ak$k\_mappedspot.L
grep -P "\t256\tMN00185:" $sam | cut -f 1,2,3 >> $seq\_L2R1_ak$k\_mappedspot.L
n1=`cat $seq\_L2R1_ak$k\_mappedspot.L | wc -l`
hifia_1 $seq\_L2R1_ak$k\_mappedspot.L $n1 > $seq\_L2R1_ak$k\_hifia_1.o
n2=`cat $seq\_L2R1_ak$k\_hifia_1.o | wc -l`
echo $seq $k $n1 $n2
done
###
flowcell=AAAHT5KHV
k=70
# i: two lanes, lane 1 and lane 2
# j: two surface, top and bottm
for i in 1 2
do
for j in 1 2
do
seq=L1R1Uniq_$i$j
n1=`cat $seq\_L2R1_ak$k\_mappedspot.L | wc -l`
n2=`cat $seq\_L2R1_ak$k\_hifia_1.o | wc -l`
if [ -e hifia_2_$seq\_L2R1_ak$k\_mappedspot_out ]
then
rm hifia_2_$seq\_L2R1_ak$k\_mappedspot_out
fi
for column in 1 2 3 4 5 6
do
for row in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do
tile=$flowcell:$i:$j$column$row
hifia_2 $seq\_L2R1_ak$k\_mappedspot.L $n1 $seq\_L2R1_ak$k\_hifia_1.o $n2 $tile >> hifia_2_$seq\_L2R1_ak$k\_mappedspot_out 2>>anye
done
done
date
done
done
####################################################

###################################################################
###################################################################
# Read 2
#
prefix=L2R2_010_

hg38=$mwd/imc/HG38
STAR --runThreadN 32 --genomeDir $hg38 --readFilesIn $L2R2 --quantMode GeneCounts --readFilesCommand zcat \ 
--outFileNamePrefix $prefix \ 
--outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1 > starlogo 2>starlogoe
 
date
grep "SN:" $prefix\Aligned.out.sam > $prefix\Aligned.NH1.sam
grep ":STAR" $prefix\Aligned.out.sam >> $prefix\Aligned.NH1.sam 
grep -P "NH:i:1\t" $prefix\Aligned.out.sam >> $prefix\Aligned.NH1.sam 
date

filetag=$prefix\Aligned.NH1
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
gtf=$mwd/genome/release104/Homo_sapiens.GRCh38.104.gtf
getgenefromgtf.pl $gtf ENSG > genensg104.b 2>>anye
cat genensg104.b | perl -p -e "s/:/\t/g" | cut -f 1,2,3,4 > genensg104clean.b
bedtools intersect -a $bam -b genensg104clean.b -wb -bed > $genicreadfile
cut -f 1,2,3,4,16 $genicreadfile > $filetag\gene5columns.L


#################################################################

n1=`cat L1R1Uniq_11_L2R1_ak70_mappedspot.L | wc -l`

# for mouse, it is ENSMUSG00000099364

grep -v "ENSG00000230876" L2R2_010_Aligned.NH1gene5columns.L > G
n2=`cat G | wc -l`

# g=ENSMUSG00000023886
# grep $g L2R2_015_Aligned.NH1gene5columns.L > hifi2$g.L
# n2=`cat hifi2$g.L | wc -l`

for i in  1 2 3 4 5 6;
do
for j in 01 02 03 04 05 06 07 08 09 10 11 12 13 14;
do
nohup ./hifiapg L1R1Uniq_11_L2R1_ak70_mappedspot.L $n1 AAAHT3CHV:1:1$i$j G $n2 > genicspotTile1$i$j.o 2>>anye &
done
done

cat genicspotTile1* > genicspotovertiles.o


####################################################
####################################################



seq=L1R1Uniq_11

for k in 40 55 70
do
bwa mem L2R1 ../lib1/raw4/$seq.fasta -a -k $k -t 64 > $seq\_L2R1_ak$k.sam 2>$seq\_L2R1_ak$k.same
sam=$seq\_L2R1_ak$k.sam

grep -P "\t0\tMN00185:" $sam | cut -f 1,2,3 > $seq\_L2R1_ak$k\_mappedspot.L
grep -P "\t256\tMN00185:" $sam | cut -f 1,2,3 >> $seq\_L2R1_ak$k\_mappedspot.L

n1=`cat $seq\_L2R1_ak$k\_mappedspot.L | wc -l`

hifia_1 $seq\_L2R1_ak$k\_mappedspot.L $n1 > $seq\_L2R1_ak$k\_hifia_1.o

n2=`cat $seq\_L2R1_ak$k\_hifia_1.o | wc -l`
echo $seq $k $n1 $n2
done

# g++ -o hifia_1 HiFianalysis_nspot_per_hifi.cpp
# g++ -o hifia_2 HiFianalysis_nhifi_per_tile.cpp
#################################################################

k=70
L=2
seq=L1R1Uniq_$L\1
n1=`cat $seq\_L2R1_ak$k\_mappedspot.L | wc -l`
n2=`cat $seq\_L2R1_ak$k\_hifia_1.o | wc -l`

for i in 1 2 3 4 5 6; 
do 
for j in 01 02 03 04 05 06 07 08 09 10 11 12 13 14; 
do 
hifia_2 $seq\_L2R1_ak$k\_mappedspot.L $n1 $seq\_L2R1_ak$k\_hifia_1.o $n2 AAAHT3CHV:$L:1$i$j >> hifiOn$seq\_k$k.L; 
done; 
date; 
done
# the output in hifiOn$seq\_k$k.L
# the sum of evenly-distributed HiFi read count per tile
# the num of data rows per tile
# the num of uniquely-resolved hifi reads per tile
# the num of all resolved hifi per tile
# the num of all mapped barcode per tile
# cout << surf << "\t" << x2 << "\t" << nline << "\t" << uniqmappedhifi2tile.size() << "\t" << hifi2tile.size() << "\t" << spot2tile.size() << "\n";

#################################################################
# /mnt/extraids/OceanStor-0/linpei/hifi/data_12/lib2
mappedbarcode_per_tile.R
# Output: hifiOnL1R1Uniq_21_k70.pdf

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

L=1
i=1
for j in 1 2 3 4 5 6; 
do 
for y in 01 02 03 04 05 06 07 08 09 10 11 12 13 14; 
do 
tile=AAAHT3CHV:$L:$i$j$y; 
nohup ./hifia_2 $seq\_L2R1_ak$k\_mappedspot.L $n1 $seq\_L2R1_ak$k\_hifia_1.o $n2 $tile > hifia_2_1_$i$j$y.o 2>>hifia_2e & 
done 
done
###################################################################
###################################################################
# Read 2
#
hg38=$mwd/imc/HG38
STAR --runThreadN 32 --genomeDir $hg38 --readFilesIn $L2R2 --quantMode GeneCounts --readFilesCommand zcat \ 
--outFileNamePrefix L2R2_010_ \ 
--outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1 > starlogo 2>starlogoe
 
date
grep "SN:" L2R2_010_Aligned.out.sam > L2R2_010_Aligned.NH1.sam
grep ":STAR" L2R2_010_Aligned.out.sam >> L2R2_010_Aligned.NH1.sam 
grep -P "NH:i:1\t" L2R2_010_Aligned.out.sam >> L2R2_010_Aligned.NH1.sam 
date

filetag=L2R2_010_Aligned.NH1
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
gtf=$mwd/genome/release104/Homo_sapiens.GRCh38.104.gtf
getgenefromgtf.pl $gtf ENSG > genensg104.b 2>>anye
cat genensg104.b | perl -p -e "s/:/\t/g" | cut -f 1,2,3,4 > genensg104clean.b
bedtools intersect -a $bam -b genensg104clean.b -wb -bed > $genicreadfile
cut -f 1,2,3,4,16 $filetag\gene.L > $filetag\gene5columns.L


#################################################################



# cat job6.sh
L=$1
i=$2
n1=$3
n2=$4

tile=AAAHT3CHV:$L:$i:
k=50
./hifia_2 L1R1UNIQ_L2R1_ak$k\_mappedspot.L $n1 L1R1UNIQ_L2R1_ak$k\_nspot_per_hifi.o $n2 $tile

## ./hifia_2 L1R1UNIQ_L2R1_ak$k\_mappedspot.L $n L1R1UNIQ_L2R1_ak$k\_nspot_per_hifi.o $n2 AAAHT3CHV:1:1114 2>some



L=$1
i=$2
n=$3
tile=AAAHT3CHV:$L:$i:
# n=`cat L1R1UNIQ_L2R1_ak50_mappedspot.L | wc -l`
date;
./hifia1 L1R1UNIQ_L2R1_ak70_mappedspot.L $n $tile > output_$L\_$i\.o
date





date
### Read 1
for i in `grep "L001" ../lib1/raw4/filelistL1R1.L`;
do
L1R1=../lib1/raw4/$i\_R1_001.fastq.gz
# i=CZ929_S3_L001
k=90
date;
bwa mem L2R1 $L1R1 -a -k $k -t 48 > $i\_L2R1_ak$k.sam 2>$i\_L2R1_ak$k.same;
date
done

if [ -e spotoL2R1_ak90_0256.cleansam ]
then
rm spotoL2R1_ak90_0256.cleansam
fi

# Tue May 31 14:45:24 PDT 2022 roughly time spending
for j in `ls *L2R1_ak90.sam`;
do
grep -P "\t0\tMN00185:" $j | cut -f 1,2,3 >> spotoL2R1_ak90_0256.cleansam
grep -P "\t256\tMN00185:" $j | cut -f 1,2,3 >> spotoL2R1_ak90_0256.cleansam
done
# Tue May 31 15:21:40 PDT 2022 (roughly time spending)

n=`cat spotoL2R1_ak90_0256.cleansam | wc -l`
./hifia1 spotoL2R1_ak90_0256.cleansam $n


#for i in `cut -f 1 ../lib1/raw4/filelistL1R1.L`; 
#do 
#k=90; 
#grep -P "\t0\tMN00185:" $i\_L2R1_ak$k.sam | cut -f 1,2,3 > spotoL2R1_ak90_0256$i.cleansam; 
#grep -P "\t256\tMN00185:" $i\_L2R1_ak$k.sam | cut -f 1,2,3 >> spotoL2R1_ak90_0256$i.cleansam; 
#n=`cat spotoL2R1_ak90_0256$i.cleansam | wc -l`; 
#nohup ./hifia1 spotoL2R1_ak90_0256$i.cleansam $n > hifia1_$i.o 2>hifia1_$i.e &
#echo $i
#date
#done



# /mnt/extraids/OceanStor-0/linpei/hifi/data_8/lib2
L=1
date
grep -P "\t0\tVH" L2R1bwatoL1R1LANE$L\_ak80.sam | cut -f 1,2,3 > L2R1bwatoL1R1LANE$L\_ak80_0256.cleansam
grep -P "\t256\tVH" L2R1bwatoL1R1LANE$L\_ak80.sam | cut -f 1,2,3 >> L2R1bwatoL1R1LANE$L\_ak80_0256.cleansam
date
##################################################################
# align HiFi R2 reads to genome
L2R2=/mnt/extraids/OceanStor-0/linpei/hifi/data_8/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz
bwa mem $mwd/genome/release105/DNAMM39 $L2R2 -a -t 64 > bwaL2R2tomm39.sam 2>>anye

filetag=bwaL2R2tomm39
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
# getgenefromgtf.pl Mus_musculus.GRCm39.105.gtf ENSMUSG > genensmusg105.b 2>>anye
cat genensmusg105.b | perl -p -e "s/:/\t/g" | cut -f 1,2,3,4 > genensmusg105clean.b
bedtools intersect -a $bam -b genensmusg105clean.b -wb -bed > $genicreadfile
cut -f 1,2,3,4,16 bwaL2R2tomm39gene.L > bwaL2R2tomm39gene.fivecolumn.L

##################################################################
file1=L2R1bwatoL1R1LANE1_ak80_0256.cleansam 
n1=`cat $file1 | wc -l`
file2=bwaL2R2tomm39gene.fivecolumn.L
n2=`cat $file2 | wc -l`
file3=$mwd/genome/release105/pericytes_cluster_genenamelist_2.L
n3=`cat $file3 | wc -l`
##################################################################

# g++ -o hifia HiFianalysis3.cpp -lz
# linpei@sysbiocomp:/mnt/extraids/OceanStor-0/linpei/hifi/data_8/lib2$ date;
./hifia L2R1bwatoL1R1LANE1_ak80_0256.cleansam $n1 bwaL2R2tomm39gene.fivecolumn.L $n2 > tmpo;date
# Sat May 21 11:00:40 PDT 2022
# Sat May 21 11:06:27 PDT 2022
date
cut -f 1 L2R1bwatoL1R1LANE$L\_ak80_0256.cleansam | sort | uniq -c > L2R1bwatoL1R1LANE$L\_ak80_0256_hifi2spot.o


# cp /mnt/extraids/OceanStor-0/linpei/hifi/data_10/HiFianalysis2.cpp HiFianalysis3.cpp






# For a HiFi library of 17,523,315 read pairs (R1 length: 100 bps)
# it takes bwa 35 mins to run the indexing process. 
#
# it took ~80 mins to finish mapping Library 1 R1 reads (NextSeq 2000, Lane 1) to our Library 2 R1 by bwa-mem using 48 thread.
#
# it took ~8h to read and handle all the raw Library 1 R1 reads by Perl and ~2h by c++
# (1,427,451,711 raw reads,101 bps long, NextSeq 2000, two lanes) by 1 thread.
# after collapased redundant Library 1 R1 reads, a total of 505,137,607 non-redundant L1R1 reads were found.
# among these, 337,392,125 L1R1 reads occured only once in the raw sequencing data.

L2R1=../raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
bwa index -p L2R1 $L2R1 > bwaindexo 2>bwaindexe

# linpei@sysbiocomp:/mnt/extraids/OceanStor-0/linpei/hifi/data_8/lib2/raw3$ grep "|1|L" wholefastq_R1_unique.L > wholefastq_R1nonredundant.L
# linpei@sysbiocomp:/mnt/extraids/OceanStor-0/linpei/hifi/data_8/lib2/raw3$ ls -cltr wholefastq_R1_unique.L wholefastq_R1nonredundant.L
# -rw-rw-r-- 1 linpei linpei 100980703760 Mar 24 14:05 wholefastq_R1_unique.L
# -rw-rw-r-- 1 linpei linpei  20800191594 Mar 24  2022 wholefastq_R1nonredundant.L





L=$1
k=$2
################################################################## the frst step: align HiFi Reads to Library 1 raw reads
for j in `cut -f 1 ur.L`
do
date
L1R1=../../lib1/$j\_L00$L\_R1_001.fastq.gz
# -a: report all alignments. this is necessary as multiple L2 reads (R1) could be aligned at one coordinate, whose L1 read used as query for BWA.
# -k: Minimum seed length. Matches shorter than INT will be missed.default: 19 bps.
bwa mem L2R1 $L1R1 -a -k $k -t 48 > $j\_L00$L\_R1_ak$k.sam 2>$j\_L00$L\_R1_ak$k\bwae
echo "Done" $j $L
done
##################################################################
# count the number of HiFi read pairs per tile:

i=$1
k=$2
echo "MN00185" > hifiread_on_Tile$i.L
for j in `cut -f 1 ur.L`
do
grep -P "AAAHT5KHV:1:$i:" ../raw4/$j\_L001_R1_ak$k\.sam | grep -P "\t0\tMN00185" | cut -f 3 >> hifiread_on_Tile$i.L
grep -P "AAAHT5KHV:1:$i:" ../raw4/$j\_L001_R1_ak$k\.sam | grep -P "\t256\tMN00185" | cut -f 3 >> hifiread_on_Tile$i.L
done

sort hifiread_on_Tile$i.L | uniq > hifiread_on_Tile$i\_uniq.L
rm hifiread_on_Tile$i.L

# count the number of HiFi read pairs on one lane:

k=80
echo "MN00185" > hifiread_on_Lane.L

for j in `cut -f 1 ur.L`
do
  grep -P "\t0\tMN00185:" ../raw4/$j\_L001_R1_ak$k\.sam | cut -f 3 >> hifiread_on_Lane.L
grep -P "\t256\tMN00185:" ../raw4/$j\_L001_R1_ak$k\.sam | cut -f 3 >> hifiread_on_Lane.L
done

sort hifiread_on_Lane.L | uniq > hifiread_on_Lane_uniq.L
  rm hifiread_on_Lane.L

# count the number of HiFi-resolved spots on tile:
# this can be merged with codes above into one shell script.
# by definition, we only consider lines with bit score = 0

i=$1
k=$2

if [ -e hifispot_on_Tile$i.L ]
then
rm hifispot_on_Tile$i.L
fi

for j in `cut -f 1 ur.L`
do
grep -P "AAAHT5KHV:1:$i:" ../raw4/$j\_L001_R1_ak$k\.sam | grep -P "\t0\tMN00185" | cut -f 1 >> hifispot_on_Tile$i.L
done

sort hifispot_on_Tile$i.L | uniq > hifispot_on_Tile$i\_uniq.L
  rm hifispot_on_Tile$i.L
  
cut -f 1 hifispot_on_Tile$i\_uniq.L | perl -p -e "s/:/\t/g" | cut -f 6,7 > hifispot_on_Tile$i\_uniq.xy
n=`cat hifispot_on_Tile$i\_uniq.L | wc -l`
calcdist hifispot_on_Tile$i\_uniq.xy $n 0 $n > hifispot_on_Tile$i\_uniq.dist;
rm hifispot_on_Tile$i\_uniq.xy

##################################################################
# align HiFi R2 reads to genome
L2R2=/mnt/extraids/OceanStor-0/linpei/hifi/data_8/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz
bwa mem $mwd/genome/release105/DNAMM39 $L2R2 -a -t 64 > bwaL2R2tomm39.sam 2>>anye

filetag=bwaL2R2tomm39
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
# getgenefromgtf.pl Mus_musculus.GRCm39.105.gtf ENSMUSG > genensmusg105.b 2>>anye
cat genensmusg105.b | perl -p -e "s/:/\t/g" | cut -f 1,2,3,4 > genensmusg105clean.b
bedtools intersect -a $bam -b genensmusg105clean.b -wb -bed > $genicreadfile
cut -f 1,2,3,4,16 bwaL2R2tomm39gene.L > bwaL2R2tomm39gene.fivecolumn.L

##################################################################











##################################################################
#####
####
###
##
#
# the script, getBWAL2R1uniqcoord.pl, uses the second column (flag) in sam for information.
# https://en.wikipedia.org/wiki/SAM_(file_format)#Bitwise_flags
# if the value in second column is 0 or 256, the L1 reads in the first column would be assigned to the L2 reads in the third column.
# a L2 read would be outputed if it has only 1 coordinates.
#
# for j in `cut -f 1 ur.L`
# do
# nohup getBWAL2R1uniqcoord.pl $j\_L00$L\_R1_ak$k.sam > $j\_L00$L\_L2R1bycoord_ak$k\o 2>>anye &
# done
#
# j=Undetermined_S0
# nohup getBWAL2R1uniqcoord.pl $j\_L00$L\_R1_ak$k.sam > $j\_L00$L\_L2R1bycoord_ak$k\o 2>>anye &
#
# n=`ps x | grep "getBWAL2R1uniqcoord.pl" | wc -l`
# echo $n
# while [ $n -gt 1 ]
# do
# sleep 600
# date
# n=`ps x | grep "getBWAL2R1uniqcoord.pl" | wc -l`
# done
#
# cat Undetermined_S0_L00$L\_L2R1bycoord_ak$k\o > L2R1bycoordL00$L\_ak$k\o
# for j in `cut -f 1 ur.L`;
# do
# cat $j\_L00$L\_L2R1bycoord_ak$k\o >> L2R1bycoordL00$L\_ak$k\o 
# done
#
#
# the below cmd run quickly. 
# cut -f 1 L2R1bycoordL001_ak95o > L2R1bycoordL001_ak95_column1_o
# sort L2R1bycoordL001_ak95_column1_o > L2R1bycoordL001_ak95_column1_sortedo
# uniq -c L2R1bycoordL001_ak95_column1_sortedo > L2R1bycoordL001_ak95_column1_sorted_uniqized.o
#
# grep -P "^\s+1\sMN00185:" L2R1bycoordL001_ak95_column1_sorted_uniqized.o | wc -l
# 1382471
#

# old script generate a large file which record a large sparse matrix 
# nohup testGeneToTileByL2_step0.pl L2R1bycoordL00$L\_ak$k\o $genicreadfile > MatrixL00$L\_geneXtile.o 2>MatrixL00$L\_geneXtile.errmsg &
#
# input:
# L2R1bycoordL00$L\_ak$k\o: a mapping from hifi reads to L1 coordinates
# $genicreadfile: a mapping from hifi reads to genes
#
# 2 output files: MatrixL00$L\_L2ReadXtile.mtx MatrixL00$L\_geneXtile.mtx
# testGeneToTile_writemtx.pl L2R1bycoordL00$L\_ak$k\o $genicreadfile MatrixL00$L\_L2ReadXtile.mtx MatrixL00$L\_geneXtile.mtx
#
# testGeneToTile_writemtx_twolanes.pl L2R1bycoordL001_ak80o L2R1bycoordL002_ak80o $genicreadfile L2ReadXtile.mtx L2ReadXtile_suppl.mtx 
#
##### parameter: tile
#### testGeneToTile_writemtx_toplot will get coord on the given tile and give color to the top-10 genes which had highest number of spots on the tile
###
##
#
# testGeneToTile_writemtx_toplot.py MatrixL00$L\_L2ReadXtile.mtx T1208 > gene2spot 2>top10gene_for_legend
# testGeneToTile_writemtx_toplot.py MatrixL00$L\_L2ReadXtile.mtx T1208 A_FILE_TO_MAP_ENSEMBL_GENE_ID_TO_GENE_SYMBOL > tmpo 2>tmpe
#  testGeneToTile_writemtx_toplot.py L2ReadXtile.mtx L1 T$tile mousegenename.L > L2ReadonLane1T$tile\o 2>L2ReadonLane1T$tile\e &
#
# testGeneToTile_plot_genes_on_1_tile.R
#
##
###
# n1=`ps x | grep "testGeneToTileByL2_step0.pl" | wc -l`
# while [ $n1 -gt 1 ]
# do
# sleep 120
# date
# n1=`ps x | grep "testGeneToTileByL2_step0.pl" | wc -l`
# done
#
# for ng in 500 1000 2000;
# do
# filein=MatrixL00$L\_L2ReadXtile.mtx
# nohup R CMD BATCH --no-save --no-restore "--args $filein 8 $ng 64" testGeneToTile_1.R 2>>testGeneToTile_1.Rerr &
# done