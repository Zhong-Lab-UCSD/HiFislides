# do the de-duplication of raw reads from recycled flow cell.
# then the unique raw reads from recycled flow cell would be used as spatial barcodes.

# working dir:
#

flowcell=AAAL33WM5

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


# working dir:
#

L2R1=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
L2R2=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz

flowcell=AAAL33WM5
seq=L1R1Uniq_11
k=40

bwa index -p L2R1 $L2R1 > bwaindexL2R1o 2>bwaindexL2R1e


k=40
seq=L1R1Uniq_11
sam=$seq\_L2R1_ak$k.sam
bwa mem L2R1 $seq.fasta -a -k $k -t 64 > $sam 2>$seq\_L2R1_ak$k.same

grep -P "\t0\tMN00185:" $sam | cut -f 1,2,3 > $seq\_L2R1_ak$k\_mappedspot.L
grep -P "\t256\tMN00185:" $sam | cut -f 1,2,3 >> $seq\_L2R1_ak$k\_mappedspot.L
