L=$1
k=$2
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

R CMD BATCH --no-save --no-restore "--args MatrixL00$L\_geneXtile.o" test2.R
