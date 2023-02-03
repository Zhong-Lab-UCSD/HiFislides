chmod u+x /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/bin/hifi_wrapper.sh
conda activate hifi

# /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/MiniSeq/data_26_13aughifi

/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/bin/hifi_wrapper.sh \
-b /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/bin \
-i /dataOS/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/STARindex_withSJ \
-g /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gtf \
-N data_26_13aughifi_wrapper \
-S /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAAL235M5 \
-1 /mnt/extraids/SDSC_NFS/linpei/hifi/data_26_13aughifi/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz \
-2 /mnt/extraids/SDSC_NFS/linpei/hifi/data_26_13aughifi/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz \
-t 32 \
-o /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/MiniSeq
