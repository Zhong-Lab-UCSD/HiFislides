# HiFi-Slides Spatial Sequencing
Contributions by Stanley, Ekko, Pei, Riccardo, Brianne in Zhong Lab.

# About  
HiFi-Slide is a super-resolution spatial transcriptomics sequencing technology. This technique captures and spatially resolves genome-wide RNA expression in a submicron resolution for fresh-frozen tissue.  

# Computational tools
All published computational tools used for HiFi-Slide data analysis are free software. We installed all of them under one conda environment.   

**Required software:**  
- [PEAR - Paired-End reAd mergeR](https://cme.h-its.org/exelixis/web/software/pear/)  
- [fastp](https://github.com/OpenGene/fastp)  
- [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)  
- [STAR](https://github.com/alexdobin/STAR)  
- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)  

**R packages:**  
- [ComplexHeatmap](http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)  

# References
- PEAR: a fast and accurate Illumina Paired-End reAd mergeR
Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593  
- Rozowsky, J., Kitchen, R. R., Park, J. J., Galeev, T. R., Diao, J., Warrell, J., Thistlethwaite, W., Subramanian, S. L., Milosavljevic, A., & Gerstein, M. (2019). exceRpt: A Comprehensive Analytic Platform for Extracellular RNA Profiling. Cell systems, 8(4), 352–357.e3. https://doi.org/10.1016/j.cels.2019.03.004

# Workflow

## 1. Deduplication of spatial barcodes (L1R1)
```
surfdedup AAAL33WM5:1:1 *_L1R1.fastq.gz > L1R1_dedup.fasta 2>L1R1_dup.txt
```
**Arguments**  
1. Identifier of the surface on recycled flowcell. For example, AAAL33WM5:1:1 indicates that the deduplication will be performed on the top surface (surface 1) of lane 1 of flowcell AAAL33WM5. 
2. The filename identifier of several L1R1 fastq files.

**Purpose**  
Remove redundant spatial barcodes. This script reads raw reads from recycled flow cell and adds "_N" to the identifier of each read, where N is the number of times the read is found on the surface (i.e. one spatial barcode could be found at N different coordinates on the surface). This script considered only reads from the surface identified by argument 1. 

**Output**  
The output of surfdedup includes two files:

1. `L1R1_dedup.fasta`: fasta file with deduplicated read sequences (obtained from STDOUT). 
2. `L1R1_dup.txt`: txt file listing all the read identifiers that shared the same read sequence. Note that when N reads shared the same sequence, only 1 of the N read identifiers would be randomly chosen and printed to `L1R1_dedup.fasta` while the remaining N - 1 read identifiers would be shown in N - 1 rows in `L1R1_dup.txt`.
  

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
hifislida.pl output_sam_file_by_BWA
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
hifislida2.pl Ouput_from_hifislida output_sam_file_by_BWA
```

**3-1. Arguments**  
argument \#1:  Output file produced by hifislida.pl
argument \#2:  Outout file produced by BWA

**3-2. Purpose**  
Count the number of HIFISLIDE reads per tile. A total of 6 X 11 tiles were available on Nextseq flowcell (sometimes it could be 6 X 14). We hypothesized that spatial barcodes on tiles coverred by tissue should be mapped with more HIFISLIDE R1 thans spatial barcodes outside tissue cover region. To this end, we count the number of HIFISLIDE R1 reads per tile. To find a simplilified solution, we only considered HIFISLIDE R1 which had only one unique spatial barcode with highest alignment score on the surface.       
**3-3. Output format**  
Column 1 - Tile ID.  
Column 2 - Number of spatially resolved HIFISLIDE R1 reads per tile (ranked in descending order)

If tiles with highest number of HIFISLIDE R1 reads tend to be located in proximity, that may indicate these tiles were covered by the tissue.

## 4. preprocessing of HIFISLIDE R2 reads  
By design, HIFISLIDE R2 sequenced the tissue RNA. It is the RNA end. In practice, one issue was the read throught by HIFISLIDE R2 into the spatial barcode. If occurred, HIFISLIDE R2 could carry sequence of the R1 from the recycled flowcell. To identify such cases, we search for the illumina R1 primer in HIFISLIDE R2 and also search for the overlap between HIFISLIDE R1 and R2 per read pair. The latter task was performed by PEAR v0.9.6 using default parameters. We excluded HIFISLIDE R2 that overlap with HIFISLIDE R1 or mapped with illumina R1 primer.  
Next, we used ``fastp`` to further process HIFISLIDE R2 reads. We had 6 different options for processing.


Option 1
```
fastp -i L2R2_1x2.fastq -o L2R2_1x2_trim_tail1_1.fastq --trim_tail1 80 --disable_quality_filtering --thread 16 > someo 2>somee;date
```

Option 2
```
fastp -i L2R2_1x2.fastq -o L2R2_1x2_trim_front1_1.fastq --trim_front1 60 --disable_quality_filtering --thread 16 > someo 2>somee;date
```

Option 3
```
p=L2R2_1x2_processed_Q0
fastp -i L2R2_1x2.fastq -o $p.fastq --disable_quality_filtering --trim_poly_g --trim_poly_x --thread 16 > $p\o 2>$p\e
```

Option 4
```
p=L2R2_1x2_processed_Q1
fastp -i L2R2_1x2.fastq -o $p.fastq --disable_quality_filtering --trim_poly_g --trim_poly_x --cut_front --cut_tail --thread 16 > $p\o 2>$p\e
```

Option 5
```
p=L2R2_1x2_processed_Q2
fastp -i L2R2_1x2.fastq -o $p.fastq --disable_quality_filtering --trim_poly_g --trim_poly_x --cut_front --thread 16 > $p\o 2>$p\e
```

Option 6
```
p=L2R2_1x2_processed_Q3
fastp -i L2R2_1x2.fastq -o $p.fastq --disable_quality_filtering --trim_poly_g --trim_poly_x --cut_tail --thread 16 > $p\o 2>$p\e
```


Here ``L2R2_1x2.fastq`` is the fastq of filtered HIFISLIDE R2 reads that not overlapped with HIFISLIDE R1 and not mapped with illumina R1 reads primer. Processed reads were then mapped to human genome using STAR or mapped to human transcriptome using BOWTIE2.  


## 5. annotate HIFISLIDE R2 reads by genes

Two different strategies were applied.  
(1) we used STAR to align HIFISLIDE R2 reads to genome and then used bedtools to obtain annotated genes per HIFISLIDE-mapped genomic locus.  
(2) we used BOWTIE2 directly map HIFISLIDE R2 to transcriptome.
For STAR usage, we set --outFilterScoreMinOverLread and --outFilterMatchNminOverLread to be 0 as [SeqScope](https://github.com/leeju-umich/Cho_Xi_Seqscope/blob/main/script/align.sh).
For BOWTIE2, we used default setting with the local alignment mode.  
If a HIFISLIDE R2 read could be mapped to a gene using one or both strategies,it would be counted for that gene.


## 6. Integrate spatial coordinates and gene information for each HiFi read pairs.
Results from Step 3 and 5 would be intergrated to provide gene annotation for each spatially resolved HIFISLIDE read pair. The output of this step would be   
column 1 - Tile ID  
column 2 - X-coord  
column 3 - Y-coord  
column 4 - HIFISLIDE Read ID  
column 5 - Gene  