# Part 1
# do the de-duplication of raw reads from recycled flow cell.
# then the unique raw reads from recycled flow cell would be used as spatial barcodes.

# working dir:
# /mnt/extraids/OceanStor-0/linpei/hifi/data_12/lib1/raw2 

flowcell=AAAL33WM5

# i=1 means Lane 1
# j=1 means top surface. (1:top, 2:bottom)
# both surfdedup and finduniqread.pl are executable programs at /home/linpei/bin

i=1
j=1
surf=$flowcell:$i:$j
surfdedup $surf *_L00$i\_R1_001.fastq.gz > L1R1Dedup_$i$j.fasta
finduniqread.pl L1R1Dedup_$i$j.fasta > /mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib1/L1R1Uniq_$i$j.fasta

# Part 2
# working dir:
# /mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2

L2R1=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
L2R2=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz

flowcell=AAAL33WM5
seq=L1R1Uniq_11
k=40

### work with the coordinate end, i.e., R1 of HiFi slides read pairs
bwa index -p L2R1 $L2R1 > bwaindexL2R1o 2>bwaindexL2R1e

i=1
j=1
spatialbarcode=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/L1R1Uniq_$i$j.fasta

sam=$seq\_L2R1_ak$k.sam
bwa mem L2R1 $spatialbarcode -a -k $k -t 64 > $sam 2>$seq\_L2R1_ak$k.same

grep -P "\t0\tMN00185:" $sam | cut -f 1,2,3 > $seq\_L2R1_ak$k\_mappedspot.L
grep -P "\t256\tMN00185:" $sam | cut -f 1,2,3 >> $seq\_L2R1_ak$k\_mappedspot.L

### work with the RNA end, i.e., R2 of HiFi slides read pairs
# mwd=/mnt/extraids/OceanStor-0/linpei
hg38=$mwd/imc/HG38
prefix=L2R2_000_
STAR --runThreadN 32 --genomeDir $hg38 --readFilesIn $L2R2 --quantMode GeneCounts --readFilesCommand zcat \
--outFileNamePrefix $prefix \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 > starlogo 2>starlogoe

# Reference:
# https://github.com/alexdobin/STAR/issues/415
# https://github.com/leeju-umich/Cho_Xi_Seqscope/blob/main/script/align.sh 

### 

