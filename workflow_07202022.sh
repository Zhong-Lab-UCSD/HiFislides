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
