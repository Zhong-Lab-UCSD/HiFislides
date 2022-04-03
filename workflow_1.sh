# linpei@sysbiocomp:/mnt/extraids/OceanStor-0/linpei/hifi/data_8/lib2/raw3$ grep "|1|L" wholefastq_R1_unique.L > wholefastq_R1nonredundant.L
# linpei@sysbiocomp:/mnt/extraids/OceanStor-0/linpei/hifi/data_8/lib2/raw3$ ls -cltr wholefastq_R1_unique.L wholefastq_R1nonredundant.L
# -rw-rw-r-- 1 linpei linpei 100980703760 Mar 24 14:05 wholefastq_R1_unique.L
# -rw-rw-r-- 1 linpei linpei  20800191594 Mar 24  2022 wholefastq_R1nonredundant.L

L=$1
k=$2
############### the frst step: align HiFi Reads to Library 1 raw reads
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
fq1=$j\_L00$L\_R1_001.fastq
gunzip -c $fq0 > $fq1
bwa mem L2R1 $fq1 -a -k $k -t 48 > $j\_L00$L\_R1_ak$k.sam 2>$j\_L00$L\_R1_ak$k\bwae
rm $fq1
echo "Done" $j $L
#############################################
for j in `cut -f 1 ur.L`
do
nohup getBWAL2R1bytile.pl $j\_L00$L\_R1_ak$k.sam wholefastq_R1nonredundant.L > $j\_L00$L\_L2R1bytile_ak$k\o 2>>anye &
done

j=Undetermined_S0
nohup getBWAL2R1bytile.pl $j\_L00$L\_R1_ak$k.sam wholefastq_R1nonredundant.L > $j\_L00$L\_L2R1bytile_ak$k\o 2>>anye &

n=`ps x | grep "getBWAL2R1bytile.pl" | wc -l`
echo $n
while [ $n -gt 1 ]
do
sleep 600
date
n=`ps x | grep "getBWAL2R1bytile.pl" | wc -l`
done

cat Undetermined_S0_L00$L\_L2R1bytile_ak$k\o > L2R1bytileL00$L\_ak$k\o
for j in `cut -f 1 ur.L`;
do
cat $j\_L00$L\_L2R1bytile_ak$k\o >> L2R1bytileL00$L\_ak$k\o 
done

filetag=bwaL2RAW2tomm39
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
bedtools intersect -a $bam -b genensmusg105.b -wb -bed > $genicreadfile

# nohup GeneRankedBySpatialReadCount.pl L2R1bytileL00$L\_ak$k\o $genicreadfile > GenelistBySpatialRdCount.L 2>GenelistBySpatialRdCount.errmsg &
nohup testGeneToTileByL2_step0.pl L2R1bytileL00$L\_ak$k\o $genicreadfile > MatrixL00$L\_geneXtile.o 2>MatrixL00$L\_geneXtile.errmsg &
# let n=n+1

n1=`ps x | grep "testGeneToTileByL2_step0.pl" | wc -l`
while [ $n1 -gt 1 ]
do
sleep 120
date
n1=`ps x | grep "testGeneToTileByL2_step0.pl" | wc -l`
done

infile=MatrixL00$L\_geneXtile.o
# the output would be MatrixL00$L\_geneXtile.o.RData
R CMD BATCH --no-save --no-restore "--args $infile" test2.R

