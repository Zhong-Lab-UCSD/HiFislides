# **HiFi-Slides Spatial Sequencing**


Works by Stanley, Ekko, Pie, Riccardo, Brianne and many others in Zhong Lab.

# **About**  
# **Dependencies**  
# **Usage**  




# **environmental variables**
flowcell=AAAL33WM5  
i=1  
j=1  
surf=$flowcell:$i:$j  
k=40  
seq=L1R1Uniq_11  
mwd=/mnt/extraids/OceanStor-0/linpei  
hg38=$mwd/imc/HG38  
gtf=$mwd/genome/release104/Homo_sapiens.GRCh38.104.gtf  
ensgname=$mwd/genome/release104/ensg2name38104.txt  


# **custom external scripts**


**surfdedup**  
```
surfdedup AAAL33WM5:1:1 *_L00$i\_R1_001.fastq.gz  
```
(1) **Arguments**  
argument \#1: identifier of the surface. For the case AAAL33WM5:1:1, it means deduplication would be performed on the top surface of lane 1 of flowell AAAL33WM5. To use bottom surface, AAAL33WM5:1:2 need to be invoked.  
argument \#2: >= 1 fastq files of raw reads.   
(2) **Purpose**  
this script read raw reads from recycled flow cell and add "_N" to the identifier of each read. For a read whose sequence could be found N times on the surface, its identifier would be labeled with "_N". This script used only reads whose identifier contain a string matched argument \#1. 


**finduniqread.pl**  
```
finduniqread.pl Output_fasta_from_surfdedup
```
(1) **Argument**  
the output fasta from surfdedup  
(2) **Purpose**  
this script will extract and print out all raw reads whose identifiers labeled with a "_1". In this way, reads whose raw sequence found only once on the surface would be output. 

**getgenefromgtf.pl**  
create a bed format text file for gene annotation.

**hifia_asort.pl**  
```
hifia_asort.pl output_sam_file_by_BWA AAAL33WM5
```
(1) **Arguments**  
argument \#1: output sam file from BWA  
argument \#2: ID of the flowcell  
(2) **Purpose**   
This script read the sam file and output spatial barcode/coordinates that aligned to each HiFi R1 read.  
(3) **Output Format**   
column 1 - spatial barcodes  
column 2 - a useless number as placeholder  
column 3 - identifiers of HiFi read pairs  
column 4 - the number of aligned spatial coordinates for this HiFi R1 read  


**hifia_1n_marker_per_spot.pl**  
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