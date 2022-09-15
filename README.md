# **HiFi-Slides Spatial Sequencing**


Contributions by Stanley, Ekko, Pie, Riccardo, Brianne and many others in Zhong Lab.

# **About**  
# **Dependencies**  
# **Usage**  




# **environmental variables**



# **Workflow**


## deduplication of spatial barcodes
```
surfdedup AAAL33WM5:1:1 *_L00$i\_R1_001.fastq.gz  
```
**Arguments**  
argument \#1: identifier of the surface on recycled flowcell. For example, AAAL33WM5:1:1 indicates deduplication would be performed on the top surface of lane 1 of flowcell AAAL33WM5. 

argument \#2: the names of >= 1 fastq.gz files

**Purpose**  
this script read raw reads from recycled flow cell and add "_N" to the identifier of each read. For a read whose sequence could be found N times on the surface, its identifier would be labeled with "_N". This indicates that one spatial barcode could be found at N different coordinates on the surface. This script considered only reads from the surface identified by argument \#1. 

**Output**  
the output of surfdedup includes two files: (1) a fasta of deduplicated Read sequence. (2)a text file listed all read identifiers that shared the same read sequence.   
Note that when N reads shared the same sequence, only 1 of N read identifiers would be randomly chosen and printed to (1) while the remaining N - 1 read identifiers would be shown in N - 1 rows in (2)
  

## Extract gene annotation from GTF 




## Align spatial barcodes to HiFiSLIDE R1 reads

Aligner bwa was used to map HIFISLIDE R1 reads to spatial barcodes (i.e. R1 reads from recycled flowcell)

L1R1Dedup.fasta is the output fasta by surfdedup

```
bwa index L1R1 L1R1Dedup.fasta
```
HIFISLIDE_R1.fastq is the raw R1 read from HIFISLIDE sequencing.

```
bwa mem -a -k 40 -t 32 L1R1 HIFISLIDE_R1.fastq
```


```
hifia_asort.pl output_sam_file_by_BWA
```

(1) **Arguments**  
argument \#1: output sam file from BWA  


(2) **Purpose**   
This script read the sam file from bwa and count the number of spatial barcode aligned with each HIFISLIDE R1 with highest alignment score.

(3) **Output Format**   
column 1 - spatial barcodes  
column 2 - a useless number as placeholder  
column 3 - identifiers of HiFi read pairs  
column 4 - the number of aligned spatial coordinates for this HiFi R1 read  


## We integrate spatial coordinates and mapped gene for each HiFi read pairs.
```
hifia_1n_marker_per_spot.pl Output_from_hifia_asort.pl $hifi2gene $flowcell hifi2gene.G $ensgname > Output_spot_to_gene.A
```

(1)**Arguments**   
argument \#1: output from **hifia_asort.pl**   
argument \#2: a tables assign genes to HiFi read pairs.  
argument \#3: ID of the flowcell  
argument \#4: a txt provides Ensembl Gene ID for those interesting genes. In this workflow, all genes found by HIFI sequencing were used.  
argument \#5: a txt provides gene official symbl for each Ensembl Gene ID.  

(2) **Purpose**  
For each hifi read pair, **hifia_1n_marker_per_spot.pl** read the spatial barcodes aligned by R1, and genes aligned by R2. Then iterated on each spatial coordiante and print all related hifi reads and related genes, creating a spot-to-gene table.  

(3) **Output Format**  
column 1 - ID of tile  
column 2 - X axis coordinate (columns)  
column 3 - Y axis coordinate (rows)  
column 4 - ID of spatial barcode  
column 5 - Ensembl Gene ID  
column 6 - official gene symbol  
column 7 - sum of read counts that mapped to this gene at this spot  
column 8 - ID of Hifi read pair. (when multiple HiFi read pairs could be mapped to 1 gene at the same coordiante, one of them would be shown randomly)    