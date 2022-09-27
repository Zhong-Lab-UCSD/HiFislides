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

Note that here we use dummy codes just to provide and overview of the main steps of the workflow. More information can be found [here](https://docs.google.com/document/d/1MvXPgTVzzeAEnmRXDRuaJMY-U_ENorQMPzqpLVOJWA0/edit#).

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
  

## 2. Align HiFi-Slide R1 reads to deduplicated spatial barcodes

Aligner BWA was used to map HiFi-Slide R1 reads (L2R1) to deduplicated spatial barcodes `L1R1_dedup.fasta`.

First, we create bwa index from `L1R1_dedup.fasta`:
```
bwa index -p L1R1_dedup L1R1_dedup.fasta
```
Then, we align HiFi-Slide R1 reads `L2R1.fastq` to the deduplicated spatial barcodes:

```
bwa mem -a -k 40 -t 32 L1R1_dedup L2R1.fastq > L2R1__L1R1_dedup.sam 2>L2R1__L1R1_dedup.log
```

## 3. Select HiFi-Slide R1 reads with highest alignment score
```
hifislida.pl L2R1__L1R1_dedup.sam > L2R1__L1R1_dedup.hifislida.o 2>L2R1__L1R1_dedup.hifislida.e
```

**Arguments**  
1. Output SAM file from BWA.  

**Purpose**   
Read the SAM file output from BWA and collect the spatial barcodes aligned with each HiFi-Slide R1 reads at the highest alignment score. If a HiFi-Slide R1 read was aligned with N spatial barcodes that tied at the highest score, all the N spatial barcodes will be outputted except when N > 1,000. Spatial barcodes aligned with lower score would be discarded.

**Output**
Tab-separated file `L2R1__L1R1_dedup.hifislida.o` with the following columns:

- Column 1: HiFi-Slide read ID. 
- Column 2: Spatial barcode read ID. Each spatial barcode ID provide its spatial coordinates explicitly.  
- Column 3: Number of non-redundant spatial barcodes aligned at the highest score.  
- Column 4: Number of all the spatial barcodes aligned at the highest score.  
- Column 5: Highest alignment score between this HiFi-Slide R1 read and the aligned spatial barcode.


## 4. Identify the Region of Interest (ROI)

### Rank the tiles by number of HiFi-Slide read pairs

```
hifislida2.pl L2R1__L1R1_dedup.hifislida.o L2R1__L1R1_dedup.sam > L2R1__L1R1_dedup.hifislida2.o 2>L2R1__L1R1_dedup.hifislida2.e
```

**Arguments**  
1. Output file produced by hifislida.pl.
2. Output file produced by BWA.

**Purpose**  
Count the number of univocally resolved HiFi reads per tile. A total of 6 X 11 tiles are available on NextSeq flowcell (sometimes they could be 6 X 14). We hypothesize that spatial barcodes on tiles covered by tissue should be mapped with more HiFi-Slide R1 reads than spatial barcodes outside the tissue covered region. To this end, we count the number of HiFi-Slide R1 reads per tile. To have a simplilified solution, we only consider HiFi-Slide R1 reads that have only one unique spatial barcode with the highest alignment score on the surface, which means considering only rows in `L2R1__L1R1_dedup.hifislida.o` where both columns 3 and 4 are equal to 1.       

**Output**  
Tab-separated file `L2R1__L1R1_dedup.hifislida2.o` with the following columns:

- Column 1: Tile ID (with a "T" at the beginning). *Is it useful to have the T? If so, this should be always present to identify the tile ID and also in the output of hifislida3* 
- Column 2: Number of spatially resolved HiFi-Slide read pairs per tile (ranked in descending order).
- Column 3: Number of mapped barcodes on the tile used for resolving the HiFi-Slide read pairs.

### Select tiles under ROI

If tiles with high number of HiFi-Slide R1 reads tend to be located in proximity, that may indicate these tiles were covered by the tissue.

Algorithm explained in the Google Doc, script `xxx`.

**Output**
`ROI_tile_IDs.txt`: Tile IDs under ROI.


## 5. Match HiFi-Slide read pairs with spatial location 

Note: output filename to be changed to something better!!

```
hifislida3.pl L2R1__L1R1_dedup.hifislida.o ROI_tile_IDs.txt L1R1_dup.txt > hifislida3.o
```

**Arguments**  
1. Output file produced by hifislida.pl.
2. Output file produced by `xxx` with the IDs of the tiles under ROI.
3. The second output file from `surfdedup`. This file provides duplicate spatial barcodes whose IDs were not shown in the output from `hifislida.pl`.

**Purpose**  
With `hifislida.pl` we obtained aligned deduplicated spatial barcodes with the highest score for each HiFi-Slide read pair. These HiFi-Slide read pairs were considered as spatially resolved. With `hifislida2.pl` we ranked tiles by their number of spatially resolved HiFi-Slide read pairs and we could manually select a few tiles as our ROI. To integrate these results we use `hifislida3.pl` to obtain all the aligned spatial barcodes located within the ROI and print out their coordinates on each tile. Notably, when multiple spatial barcodes share the same sequence but from different coordinates, `hifislida3.pl` prints out all their coordinates.

**Output**  
Tab-separated file `hifislida3.o` with the following columns:

- Column 1: HiFi-Slide read ID.
- Column 2: Tile ID (only tiles under the ROI provided as input).
- Column 3: X-coord on the tile (columns).
- Column 4: Y-coord on the tile (rows).
- Column 5: N as the number of total spatial coordinates where this HiFi-Slide read could be aligned to spatial barcodes. This is used to "weight" HiFi-Slide reads. For example, if a HiFi-Slide read has N = 8, it would be counted as 1/8 at any of these 8 coordinates.


## 6. Preprocessing of HiFi-Slide R2 reads  

By design, HiFi-Slide R2 reads contain the sequences of the tissue RNA (RNA end of the read pair). 

It is possible that HiFi-Slide R2 sequencing reads the cDNA fragment but unintendedly reads through the fragment from one end to the other end. In other words, HiFi-Slide R2 could mistakenly include a portion of R1 in the recycled flow cell. To identify such cases, we search for overlap between L2R1 and L2R2 using the software PEAR (v0.9.6, minimum overlap size = 10 bp) and we filter out such L2R2 reads. Next, we trim the front 60 bp to remove Illumina Read 1 primers in L2R2 using the software fastp (v0.23.2). Note that besides trimming, fastp also reduces the number of reads because it filters out reads that are too short.


### Filter out overlapping L2R1 and L2R2 (pear)

```
pear \
-f L2R1 \
-r L2R2 \
-v 10 \
-j 32 \
-o L2R2_pear
```

**Arguments**  
- `-f`: Name of the file with the forward paired-end reads.
- `-r`: Name of the file with the reverse paired-end reads.
- `-v`: Minimum overlap size (default: 10).
- `-j`: Number of threads.
- `-o`: Basename of the output files.

**Output**  
The software produces several fastq files with different contents. For our purposes, we are interested in `unassembled.reverse.fastq`, which contains the L2R2 reads not overlapping L2R1.


### Trimming Illumina Read 1 primers (fastp)

The software fastp can trim only 30 bp at the time, thus it is run twice consecutively. The first input file `-i` is the output fastq from pear, the second input file is the output file from the first round of fastp. With this step we also remove reads that are too short.

```
fastp \
-i L2R2.pear_filter.fastq \
-o L2R2.trim_front_temp.fastq \
-h L2R2.trim_front_temp.log.html \
-j L2R2.trim_front_temp.log.json \
--trim_front1 30 \
--disable_quality_filtering \
--thread 16

fastp \
-i L2R2.trim_front_temp.fastq \
-o L2R2.trim_front_60.fastq \
-h L2R2.trim_front_60.log.html \
-j L2R2.trim_front_60.log.json \
--trim_front1 30 \
--disable_quality_filtering \
--thread 16
```

**Output**  
Fastq file of preprocessed L2R2 reads `L2R2.trim_front_60.fastq` that will go into the following mapping steps.


<!-- Option 1
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

Here ``L2R2_1x2.fastq`` is the fastq of filtered HiFi-Slide R2 reads that not overlapped with HiFi-Slide R1 and not mapped with illumina R1 reads primer. Processed reads were then mapped to human genome using STAR or mapped to human transcriptome using BOWTIE2.  
-->

## 7. Mapping HiFi-Slide R2 reads

1. We used STAR to align HiFi-Slide R2 reads to the genome.
2. we used BOWTIE2 to align HiFi-Slide R2 to the transcriptome.

### Mapping to the genome

For STAR, we set `--outFilterScoreMinOverLread` and `--outFilterMatchNminOverLread` to be 0 as in [SeqScope](https://github.com/leeju-umich/Cho_Xi_Seqscope/blob/main/script/align.sh).

```
STAR \
--genomeDir $STAR_INDEX \
--readFilesIn L2R2.trim_front_60.fastq \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outSAMattributes All \
--outFileNamePrefix $L2_DIR/L2R2_mapping/genome/L2R2_genome. \
--sjdbGTFfile $annotation_gtf_file \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--runThreadN 32
```

Next, uniquely mapped reads are selected using samtools.

```
samtools view -@ 32 -b -h -q 255 \
-o $L2_DIR/L2R2_mapping/genome/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
$L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam 
```

Finally, uniquely mapped reads are associated with genes using bedtools. `Homo_sapiens.GRCh38.84.chr.gene.gtf` is a modified version of the original GTF file, where only full gene body coordinates are kept (column 3 equal to "gene"), the additonal information (such as exons) were discarded. The final `HiFi_L2R2_genome.bed` shows each L2R2 read associated with the corresponding gene.

```
bedtools intersect \
-a $L2_DIR/L2R2_mapping/genome/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
-b Homo_sapiens.GRCh38.84.chr.gene.gtf \
-wb -bed | cut -f 1,2,3,4,21 > $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome.bed
```


### Mapping to the transcriptome

For BOWTIE2, we used default setting with the local alignment mode.  
If a HiFi-Slide R2 read could be mapped to a gene using one or both strategies,it would be counted for that gene.




## 8. Integrate spatial coordinates and gene information for each HiFi read pairs.
Results from Step 3 and 5 would be intergrated to provide gene annotation for each spatially resolved HiFi-Slide read pair. The output of this step would be   
column 1 - Tile ID  
column 2 - X-coord  
column 3 - Y-coord  
column 4 - HiFi-Slide Read ID  
column 5 - Gene  