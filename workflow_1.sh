
# For a HiFi library of 17,523,315 read pairs (R1 length: 100 bps)
# it takes bwa 35 mins to run the indexing process. 

# it took ~80 mins to finish mapping Library 1 R1 reads (NextSeq 2000, Lane 1) to our Library 2 R1 by bwa-mem using 48 thread.

# it took ~8h to read and handle all the raw Library 1 R1 reads (1,427,451,711 raw reads,101 bps long, NextSeq 2000, two lanes) by 1 thread.
# after collapased redundant Library 1 R1 reads, a total of 505,137,607 non-redundant L1R1 reads were found.
# among these, 337,392,125 L1R1 reads occured only once in the raw sequencing data.

L2R1RAWGZ=../raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
bwa index -p L2R1 $L2R1RAWGZ > bwaindexo 2>bwaindexe


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
fq0=../../lib1/$j\_L00$L\_R1_001.fastq.gz
# -a: report all alignments. this is necessary as multiple L2 reads (R1) could be aligned at one coordinate, whose L1 read used as query for BWA.
# -k: Minimum seed length. Matches shorter than INT will be missed.default: 19 bps.
bwa mem L2R1 $fq0 -a -k $k -t 48 > $j\_L00$L\_R1_ak$k.sam 2>$j\_L00$L\_R1_ak$k\bwae
echo "Done" $j $L
done

j=Undetermined_S0
fq0=../../lib1/$j\_L00$L\_R1_001.fastq.gz
bwa mem L2R1 $fq0 -a -k $k -t 48 > $j\_L00$L\_R1_ak$k.sam 2>$j\_L00$L\_R1_ak$k\bwae
echo "Done" $j $L

##################################################################
for j in `cut -f 1 ur.L`
do
nohup getBWAL2R1uniqcoord.pl $j\_L00$L\_R1_ak$k.sam > $j\_L00$L\_L2R1bycoord_ak$k\o 2>>anye &
done

j=Undetermined_S0
nohup getBWAL2R1uniqcoord.pl $j\_L00$L\_R1_ak$k.sam > $j\_L00$L\_L2R1bycoord_ak$k\o 2>>anye &

n=`ps x | grep "getBWAL2R1uniqcoord.pl" | wc -l`
echo $n
while [ $n -gt 1 ]
do
sleep 600
date
n=`ps x | grep "getBWAL2R1uniqcoord.pl" | wc -l`
done

cat Undetermined_S0_L00$L\_L2R1bycoord_ak$k\o > L2R1bycoordL00$L\_ak$k\o
for j in `cut -f 1 ur.L`;
do
cat $j\_L00$L\_L2R1bycoord_ak$k\o >> L2R1bycoordL00$L\_ak$k\o 
done

#
# the below cmd run quickly. 
# cut -f 1 L2R1bycoordL001_ak95o > L2R1bycoordL001_ak95_column1_o
# sort L2R1bycoordL001_ak95_column1_o > L2R1bycoordL001_ak95_column1_sortedo
# uniq -c L2R1bycoordL001_ak95_column1_sortedo > L2R1bycoordL001_ak95_column1_sorted_uniqized.o
#
# grep -P "^\s+1\sMN00185:" L2R1bycoordL001_ak95_column1_sorted_uniqized.o | wc -l
# 1382471

##################################################################
L2R2RAWGZ=../raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz
bwa mem ../../../../genome/release105/DNAMM39 $L2R2RAWGZ -t 64 > bwaL2RAW2tomm39.sam 2>>anye

filetag=bwaL2RAW2tomm39
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
# getgenefromgtf.pl Mus_musculus.GRCm39.105.gtf ENSMUSG > genensmusg105.b 2>>anye
bedtools intersect -a $bam -b genensmusg105.b -wb -bed > $genicreadfile

# old script generate a large file which record a large sparse matrix 
# nohup testGeneToTileByL2_step0.pl L2R1bycoordL00$L\_ak$k\o $genicreadfile > MatrixL00$L\_geneXtile.o 2>MatrixL00$L\_geneXtile.errmsg &
#
#
testGeneToTile_writemtx.pl L2R1bycoordL00$L\_ak$k\o $genicreadfile > MatrixL00$L\_L2ReadXtile.mtx 2>MatrixL00$L\_geneXtile.mtx

n1=`ps x | grep "testGeneToTileByL2_step0.pl" | wc -l`
while [ $n1 -gt 1 ]
do
sleep 120
date
n1=`ps x | grep "testGeneToTileByL2_step0.pl" | wc -l`
done

for ng in 500 1000 2000;
do
filein=MatrixL00$L\_geneXtile.mtx
nohup R CMD BATCH --no-save --no-restore "--args $filein 8 $ng 64" testGeneToTile_1.R 2>>testGeneToTile_1.Rerr &
done
