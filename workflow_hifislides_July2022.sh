# Part 1
# do the de-duplication of raw reads from recycled flow cell.
# then the unique raw reads from recycled flow cell would be used as spatial barcodes.

# working dir:
# /mnt/extraids/OceanStor-0/linpei/hifi/data_12/lib1/raw2 

flowcell=AAAL33WM5

# i=1 means Lane 1
# j=1 means top surface. (1:top, 2:bottom)
# both surfdedup(C++) and finduniqread.pl are executable programs at /home/linpei/bin

i=1
j=1
surf=$flowcell:$i:$j
surfdedup $surf *_L00$i\_R1_001.fastq.gz > L1R1Dedup_$i$j.fasta
finduniqread.pl L1R1Dedup_$i$j.fasta > /mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib1/L1R1Uniq_$i$j.fasta

# in L1R1Dedup_$i$j.fasta, 
# reads whose sequence was found N times on the surface were shown.
# in L1R1Uniq_$i$j.fasta, 
# only reads whose sequence was found 1 time on the surface were shown.

# Part 2
# working dir:
# /mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2
# R1 reads of HiFi sequencing - spatial end/coordiante end
L2R1=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
# R2 reads of HiFi sequencing - RNA end
L2R2=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz

### work with the coordinate end, i.e., R1 of HiFi slides read pairs
#
bwa index -p L2R1 $L2R1 > bwaindexL2R1o 2>bwaindexL2R1e

i=1
j=1
spatialbarcode=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/L1R1Uniq_$i$j.fasta

flowcell=AAAL33WM5
seq=L1R1Uniq_11
k=40
sam=$seq\_L2R1_ak$k.sam
bwa mem L2R1 $spatialbarcode -a -k $k -t 64 > $sam 2>$seq\_L2R1_ak$k.same

# these two following steps are now less useful.
grep -P "\t0\tMN00185:" $sam | cut -f 1,2,3 > $seq\_L2R1_ak$k\_mappedspot.L
grep -P "\t256\tMN00185:" $sam | cut -f 1,2,3 >> $seq\_L2R1_ak$k\_mappedspot.L

### work with the RNA end, i.e., R2 of HiFi slides read pairs
#

# 
mwd=/mnt/extraids/OceanStor-0/linpei
hg38=$mwd/imc/HG38
prefix=L2R2_000_
STAR --runThreadN 32 --genomeDir $hg38 --readFilesIn $L2R2 --quantMode GeneCounts --readFilesCommand zcat \
--outFileNamePrefix $prefix \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 > starlogo 2>starlogoe

# Reference for the usage of STAR:
# (1) https://github.com/alexdobin/STAR/issues/415
# (2) this below github page is the pipeline for seq-scope.
# https://github.com/leeju-umich/Cho_Xi_Seqscope/blob/main/script/align.sh 


grep "SN:" $prefix\Aligned.out.sam > $prefix\Aligned.NH1.sam
grep ":STAR" $prefix\Aligned.out.sam >> $prefix\Aligned.NH1.sam 
grep -P "NH:i:1\t" $prefix\Aligned.out.sam >> $prefix\Aligned.NH1.sam 

filetag=$prefix\Aligned.NH1
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
gtf=$mwd/genome/release104/Homo_sapiens.GRCh38.104.gtf

# getgenefromgtf.pl is an executable file at /home/linpei/bin
getgenefromgtf.pl $gtf ENSG > genensg104.b 2>>anye

cat genensg104.b | perl -p -e "s/:/\t/g" | cut -f 1,2,3,4 > genensg104clean.b
bedtools intersect -a $bam -b genensg104clean.b -wb -bed > $genicreadfile
cut -f 1,2,3,4,16 $genicreadfile > $filetag\gene5columns.L

# hifi2gene: assign each HiFi read pairs to genes
hifi2gene=L2R2_000_Aligned.NH1gene5columns.L

### 

# hifia_asort.pl and hifia_1n_marker_per_spot.pl are executable programs at /home/linpei/bin

seq=L1R1Uniq_11
k=40
sam=$seq\_L2R1_ak$k.sam
flowcell=AAAL33WM5

# by hifia_asort.pl, a HiFi R1 aligned to N spatial barcode (N spatial spots/coordiantes) would be counted as 1/N

# the 4 columns in $seq\_L2R1_ak$k\_mappedspot_1n.L were:
# 1) spatial barcodes
# 2) a useless number as placeholder
# 3) identifiers of HiFi read pair
# 4) the number N for this HiFi R1 read.


hifia_asort.pl $sam $flowcell > $seq\_L2R1_ak$k\_mappedspot_1n.L
cut -f 5 $hifi2gene | sort | uniq > hifi2gene.G
ensgname=/mnt/extraids/OceanStor-0/linpei/genome/release104/ensg2name38104.txt

hifia_1n_marker_per_spot.pl $seq\_L2R1_ak$k\_mappedspot_1n.L $hifi2gene $flowcell hifi2gene.G $ensgname > Output_spot_to_gene.A 
# columns in Output_spot_to_gene.A were:
(1) ID of Tile
(2) X: columns - coordinate for a spot/spatial barcode
(3) Y: rows - coordinante for a spot/spatial barcode
(4) ID of spatial barcode
(5) Ensembl Gene ID
(6) official gene symbol
(7) sum of read counts that mapped to each gene at each spot
(8) ID of Hifi read pair. Multiple HiFi read pairs could be mapped to 1 gene at the same coordiante. 
For such cases, one of these read pairs would be shown randomly.


In Output_spot_to_gene.A, 
* one coordinate(spatial barcode) could be assigned with G genes.
* only HiFi reads mapped to 1 gene were used.
if one coordinate C1 was assigned with gene G1 and G2,
there would be 2 lines about C1, one for G1 and one for G2.
if 
R2 of hifi read pairs A,B,C were mapped to G1 and 
R1 of hifi read pairs A,B,C mapped to C1,
the read count of A,B,C would be summed up and shown in column 7th.
If R1 of A,B,c were mapped to NA,NB,NC spatial coordiantes on the surface, each one would be counted as 1/NA, 1/NB and 1/NC.




