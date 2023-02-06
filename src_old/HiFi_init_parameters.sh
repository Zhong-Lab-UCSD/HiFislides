################## INPUT PARAMETERS (TO BE UPDATED PER EACH SAMPLE)
OUT_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/IGM
BIN_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/bin
N_THREADS=32

SAMPLE_NAME=HiFi_placenta_1

# Flowcell and surface identifiers
L1_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes # directory of the spatial barcodes
flowcell_type="NextSeq" # one of: MiniSeq, NextSeq
flowcell="AAALN5GM5"

# Raw reads of HiFi Slides sequencing
L2R1_FASTQ=/mnt/extraids/SDSC_NFS/rcalandrelli/Lab_sequencing/IGM/fastq/221018_A00953_0641_BHGHWNDRX2/JP_placenta_hifi_05oct22_S2_L002_R1_001.fastq.gz
L2R2_FASTQ=/mnt/extraids/SDSC_NFS/rcalandrelli/Lab_sequencing/IGM/fastq/221018_A00953_0641_BHGHWNDRX2/JP_placenta_hifi_05oct22_S2_L002_R2_001.fastq.gz

# Annotation file hg38
annotation_gtf_file=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gtf

### Genome reference indexes
STAR_INDEX=/dataOS/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/STARindex_withSJ