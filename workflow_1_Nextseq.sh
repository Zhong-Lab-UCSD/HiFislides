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

##################################################################
# align HiFi R2 reads to genome
L2R2=../raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz
bwa mem ../../../../genome/release105/DNAMM39 $L2R2 -t 64 > bwaL2RAW2tomm39.sam 2>>anye

filetag=bwaL2RAW2tomm39
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
# getgenefromgtf.pl Mus_musculus.GRCm39.105.gtf ENSMUSG > genensmusg105.b 2>>anye
bedtools intersect -a $bam -b genensmusg105.b -wb -bed > $genicreadfile

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
