flowcell=AAAL33WM5.   
i=1.   
j=1.   
surf=$flowcell:$i:$j. 

cd /mnt/extraids/OceanStor-0/linpei/hifi/data_12/lib1/raw2
surfdedup $surf *_L00$i\_R1_001.fastq.gz > L1R1Dedup_$i$j.fasta
finduniqread.pl L1R1Dedup_$i$j.fasta > /mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/L1R1Uniq_$i$j.fasta
