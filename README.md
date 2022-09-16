# **HiFi-Slides Spatial Sequencing**


Contributions by Stanley, Ekko, Pie, Riccardo, Brianne and many others in Zhong Lab.

# **About**  
HiFi-Slide sequencing is a super-resolution spatial transcriptomics sequencing technology. This technique captures and spatially resolves genome-wide RNA expression in a submicron resolution for fresh-frozen tissue.  
# **Dependencies**  
All used published tools are free softwares. We installed all of them under one conda environment.
# **References**
Rozowsky, J., Kitchen, R. R., Park, J. J., Galeev, T. R., Diao, J., Warrell, J., Thistlethwaite, W., Subramanian, S. L., Milosavljevic, A., & Gerstein, M. (2019). exceRpt: A Comprehensive Analytic Platform for Extracellular RNA Profiling. Cell systems, 8(4), 352–357.e3. https://doi.org/10.1016/j.cels.2019.03.004

# **Workflow**


## 1. deduplication of spatial barcodes
```
surfdedup AAAL33WM5:1:1 *_L00$i\_R1_001.fastq.gz  
```
**1-1. Arguments**  
argument \#1: identifier of the surface on recycled flowcell. For example, AAAL33WM5:1:1 indicates deduplication would be performed on the top surface of lane 1 of flowcell AAAL33WM5. 

argument \#2: the names of >= 1 fastq.gz files

**1-2. Purpose**  
Remove redundant spatial barcodes. This script read raw reads from recycled flow cell and add "_N" to the identifier of each read. For a read whose sequence could be found N times on the surface, its identifier would be labeled with "_N". This indicates that one spatial barcode could be found at N different coordinates on the surface. This script considered only reads from the surface identified by argument \#1. 

**1-3. Output**  
the output of surfdedup includes two files: (1) a fasta of deduplicated Read sequence. (2)a text file listed all read identifiers that shared the same read sequence.   
Note that when N reads shared the same sequence, only 1 of N read identifiers would be randomly chosen and printed to (1) while the remaining N - 1 read identifiers would be shown in N - 1 rows in (2)
  

## 2. Align HIFISLIDE R1 reads to deduplicated spatial barcodes

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

**2-1. Arguments**  
argument \#1: output sam file from BWA  

**2-2. Purpose**   
This script read the sam file output from bwa and collect spatial barcodes aligned with each HIFISLIDE R1 at highest alignment score. If a HIFISLIDE R1 was aligned with N spatial barcodes that tied at the highest score, all N spatial barcodes would be output except when N > 1,000. Spatial barcodes aligned with lower score would be discarded.

**2-3. Output Format**   
column 1 - HIFISLIDE Read ID   
column 2 - spatial barcode Read ID. Each spatial barcode ID provide its spatial coordinantes explicitly.  
column 3 - Number of non-redundant spatial barcode aligned at the highest score.  
column 4 - Number of all spatial barcode aligned at the highest score.  
column 5 - highest alignment score between this HIFISLIDE R1 and aligned spatial barcode.  

## 3. identify Region Of Interest(ROI)

```
hifislida2.pl 
```

**3-1. Arguments**  

**3-2. Purpose**  
Count the number of HIFISLIDE reads per tile. A total of 6 X 11 tiles were available on Nextseq flowcell (sometimes it could be 6 X 14). We hypothesized that spatial barcodes on tiles coverred by tissue should be mapped with more HIFISLIDE R1 thans spatial barcodes outside tissue cover region. To this end, we count the number of HIFISLIDE R1 reads per tile. To find a simplilified solution, we only considered HIFISLIDE R1 which had only one unique spatial barcode with highest alignment score on the surface.       
**3-3. Output format**  
Column 1 - Tile ID.  
Column 2 - Number of spatially resolved HIFISLIDE R1 reads per tile (ranked in descending order)

If tiles with highest number of HIFISLIDE R1 reads tend to be located in proximity, that may indicate these tiles were likely covered by the tissue.

## 4. preprocessing of HIFISLIDE R2 reads  
By design, HIFISLIDE R2 sequenced the tissue RNA. It is the RNA end. In practice, one issue was the read throught by HIFISLIDE R2 into the spatial barcode. If occurred, HIFISLIDE R2 could carry sequence of the R1 from the recycled flowcell. To identify such cases, we search for the illumina R1 primer in HIFISLIDE R2 and also search for the overlap between HIFISLIDE R1 and R2 per read pair. The latter task was performed by PEAR v0.9.6 using default parameters. We excluded HIFISLIDE R2 that overlap with HIFISLIDE R1 or mapped with illumina R1 primer.  
Next, we used ``fastp`` to further process HIFISLIDE R2 reads.

**Option 1**

```
fastp -i L2R2_1x2.fastq -o L2R2_1x2_trim_tail1_1.fastq --trim_tail1 100 --disable_quality_filtering --thread 16
```

```
fastp -i L2R2_1x2_trim_tail1_1.fastq -o L2R2_1x2_trim_tail1_2.fastq --trim_tail1 30 --disable_quality_filtering --thread 16
```
**Option 2**
```
p=L2R2_1x2_processed_Q0
fastp -i L2R2_1x2.fastq -o $p.fastq --disable_quality_filtering --trim_poly_g --trim_poly_x --thread 16 > $p\o 2>$p\e
```
**Option 3**
```
p=L2R2_1x2_processed_Q1
fastp -i L2R2_1x2.fastq -o $p.fastq --trim_poly_g --trim_poly_x --thread 16 > $p\o 2>$p\e
```
**Option 4**
```
p=L2R2_1x2_processed_Q2
fastp -i L2R2_1x2.fastq -o $p.fastq --trim_poly_g --trim_poly_x --cut_front --cut_tail --thread 16 > $p\o 2>$p\e
```
Here ``L2R2_1x2.fastq`` is the fastq of filtered HIFISLIDE R2 reads that not overlapped with HIFISLIDE R1 and not mapped with illumina R1 reads. Processed reads were then mapped to human genome using STAR or mapped to human transcriptome using BOWTIE2.  
If a HIFISLIDE R2 read could be mapped to a  gene using any of these options,that gene was assigned to the R2 read.



## 5. annotate HIFISLIDE R2 reads by genes/transcripts

Two different strategies were applied here.  
(1) we used STAR to align HIFISLIDE R2 reads to genome and then used bedtools to obtain annotated genes per HIFISLIDE-mapped genomic locus.  
(2) we used BOWTIE2 directly map HIFISLIDE R2 to transcriptome.
For STAR usage, we set --outFilterScoreMinOverLread and --outFilterMatchNminOverLread to be 0 as SeqScope.
For BOWTIE2, we used default setting with the local alignment mode.
If a HIFISLIDE R2 read could be mapped to any gene using one or both options,that gene was assigned to the R2 read.


## 6. Integrate spatial coordinates and gene information for each HiFi read pairs.