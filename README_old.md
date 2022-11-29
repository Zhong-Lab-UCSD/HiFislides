# HiFi-Slides Spatial Sequencing
Contributions by Stanley, Riccardo, "Pie", Ekko, Brianne in Zhong Lab.

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
- [Bedtools](https://bedtools.readthedocs.io/en/latest/)
- [SAMtools](https://www.htslib.org)
- [Seqtk](https://github.com/lh3/seqtk) TBD???


**R packages:**  
- [ComplexHeatmap](http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)  

# References
- PEAR: a fast and accurate Illumina Paired-End reAd mergeR
Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593  
- Rozowsky, J., Kitchen, R. R., Park, J. J., Galeev, T. R., Diao, J., Warrell, J., Thistlethwaite, W., Subramanian, S. L., Milosavljevic, A., & Gerstein, M. (2019). exceRpt: A Comprehensive Analytic Platform for Extracellular RNA Profiling. Cell systems, 8(4), 352–357.e3. https://doi.org/10.1016/j.cels.2019.03.004

# Workflow

Note that here we use dummy codes just to provide and overview of the main steps of the workflow. More information can be found [here](https://docs.google.com/document/d/1MvXPgTVzzeAEnmRXDRuaJMY-U_ENorQMPzqpLVOJWA0/edit#).

The data processing is performed separately with two shell scripts for library 1 of spatial barcodes and library 2 of HiFi-Slide read pairs:

- **`spatial_barcode_processing.sh`**: The output files from this are saved into a folder named `flowcell_ID`. The entire log of this workflow with the executed commands and the timestamps is saved into a file named `flowcell_ID.log`
- **`HiFi_processing_pipeline.sh`**: steps 1-8 below. The output files from this are saved into a folder named `sample_name`. The entire log of this workflow with the executed commands and the timestamps is saved into a file named `sample_name.log`.


## Deduplication of spatial barcodes (L1R1) and BWA index creation

This step is related to the flowcell library (library 1). It is independendent on the HiFi-Slide library.

```
flowcell_type="NextSeq"
flowcell="AAAL33WM5"

surfdedup \
$surface \
*_L1R1.fastq.gz > L1R1_dedup.fasta 2>L1R1_dup.txt
```
**Arguments**  
1. Identifier of the surface on recycled flowcell (`$surface` is automatically generated by the script using input parameters `flowcell_type` and `flowcell`). For example, `AAAL33WM5:1:1` indicates that the deduplication will be performed on the top surface (surface 1) of lane 1 of flowcell AAAL33WM5. In case of MiniSeq, we do not know if the tissue is on the top or bottom surface, therefore only `AAAL33WM5:1:` would be inputted here to consider both the surfaces.
2. The filename identifier of several L1R1 fastq files.

**Purpose**  
Remove redundant spatial barcodes. This script reads raw reads from recycled flow cell and adds "_N" to the identifier of each read, where N is the number of times the read is found on the surface (i.e. one spatial barcode could be found at N different coordinates on the surface). This script considered only reads from the surface identified by argument 1. 

**Output**  
The output of surfdedup includes two files:

1. `L1R1_dedup.fasta`: fasta file with deduplicated read sequences (obtained from STDOUT). 
2. `L1R1_dup.txt`: txt file listing all the read identifiers that shared the same read sequence (obtained from STDERR). Note that when N reads shared the same sequence, only 1 of the N read identifiers would be randomly chosen and printed to `L1R1_dedup.fasta` while the remaining N - 1 read identifiers would be shown in N - 1 rows in `L1R1_dup.txt`.

Then, we create bwa index from `L1R1_dedup.fasta`:
```
bwa index \
-p L1R1_dedup \
L1R1_dedup.fasta
```

## 1. Preprocessing of HiFi-Slide R2 reads  

By design, HiFi-Slide R2 reads contain the sequences of the tissue RNA (RNA end of the read pair). 

It is possible that HiFi-Slide R2 sequencing reads the cDNA fragment but unintendedly reads through the fragment from one end to the other end. In other words, HiFi-Slide R2 could mistakenly include a portion of R1 in the recycled flow cell. To identify such cases, we search for overlap between L2R1 and L2R2 using the software PEAR (v0.9.6, minimum overlap size = 10 bp) and we filter out such L2R2 reads. 

Next, we trim Illumina Read 1 primers, polyG and polyX tails in L2R2 using the software fastp (v0.23.2). Note that besides trimming, fastp also reduces the number of reads because it filters out reads that are too short.


### Filter out L2R2 overlapping L2R1 (PEAR)

```
pear \
-f L2R1.fastq \
-r L2R2.fastq \
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


### Trimming Illumina Read 1 primers, polyG, polyX tails (FASTP)

```
fastp \
-i L2R2.pear_filter.fastq \
-o L2R2.fastp_filter.fastq \
-h L2R2.fastp_filter.log.html \
-j L2R2.fastp_filter.log.json \
--trim_poly_g \
--trim_poly_x \
--disable_quality_filtering \
--thread 16
```

**Output**  
Fastq file of preprocessed L2R2 reads `L2R2.fastp_filter.fastq` that will go into the following mapping steps.


## 2. Mapping HiFi-Slide R2 reads

1. We used STAR to align HiFi-Slide R2 reads to the genome.
2. we used Bowtie 2 to align HiFi-Slide R2 reads to the transcriptome.

### Mapping to the genome

For STAR, we set `--outFilterScoreMinOverLread 0` and `--outFilterMatchNminOverLread 0` as in [SeqScope](https://github.com/leeju-umich/Cho_Xi_Seqscope/blob/main/script/align.sh).

```
L2R2_FILTER_FASTQ=$L2_DIR/L2R2_preprocessing/L2R2.trim_tail_60.fastq

STAR \
--genomeDir $STAR_INDEX \
--readFilesIn $L2R2_FILTER_FASTQ \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outSAMattributes All \
--outFileNamePrefix $L2_DIR/L2R2_mapping/genome/L2R2_genome. \
--sjdbGTFfile gencode.v41.annotation.gtf \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--runThreadN 32
```

Next, uniquely mapped reads are selected using samtools (NOTE: maybe here using MAPQ30 would be good enough).

```
samtools view -@ 32 -b -h -q 255 \
-o $L2_DIR/L2R2_mapping/genome/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
$L2_DIR/L2R2_mapping/genome/L2R2_genome.Aligned.sortedByCoord.out.bam 
```

Finally, uniquely mapped reads are associated with genes using bedtools. `gencode.v41.annotation.gene.gtf` is a modified version of the original GTF file `gencode.v41.annotation.gtf`, where only full gene body coordinates are kept (column 3 equal to "gene"), the additional rows (such as those related to exons for example) were discarded.

```
bedtools intersect \
-a $L2_DIR/L2R2_mapping/genome/L2R2_genome.uniquelyAligned.sortedByCoord.out.bam \
-b gencode.v41.annotation.gene.gtf \
-wb -bed | cut -f 1,2,3,4,21 > $L2_DIR/L2R2_mapping/genome/HiFi_L2R2_genome_temp.bed
```

`HiFi_L2R2_genome_temp.bed` shows each L2R2 read associated with the corresponding gene where gene information are formatted as in a GTF file. 

Example:
```
MN00185:308:000H3YMVH:1:11104:22416:20248       gene_id "ENSG00000261456.6"; gene_type "protein_coding"; gene_name "TUBB8"; level 2; hgnc_id "HGNC:20773"; havana_gene "OTTHUMG00000174803.2";
MN00185:308:000H3YMVH:1:22102:18602:1451        gene_id "ENSG00000015171.20"; gene_type "protein_coding"; gene_name "ZMYND11"; level 2; hgnc_id "HGNC:16966"; havana_gene "OTTHUMG00000017526.8";
MN00185:308:000H3YMVH:1:11103:2816:12628        gene_id "ENSG00000015171.20"; gene_type "protein_coding"; gene_name "ZMYND11"; level 2; hgnc_id "HGNC:16966"; havana_gene "OTTHUMG00000017526.8";
MN00185:308:000H3YMVH:1:12102:3409:4239         gene_id "ENSG00000015171.20"; gene_type "protein_coding"; gene_name "ZMYND11"; level 2; hgnc_id "HGNC:16966"; havana_gene "OTTHUMG00000017526.8";
```

Additional processing will transform the file into the final `HiFi_L2R2_genome.bed` with the following columns:

- Column 1: L2R2 read ID.
- Column 2: Gene ID.
- Column 3: Gene name.
- Column 4: Gene biotype.

Example:
```
MN00185:308:000H3YMVH:1:11104:22416:20248       ENSG00000261456.6       TUBB8   protein_coding
MN00185:308:000H3YMVH:1:22102:18602:1451        ENSG00000015171.20      ZMYND11 protein_coding
MN00185:308:000H3YMVH:1:11103:2816:12628        ENSG00000015171.20      ZMYND11 protein_coding
MN00185:308:000H3YMVH:1:12102:3409:4239         ENSG00000015171.20      ZMYND11 protein_coding
```


### Mapping to the transcriptome

For Bowtie 2, we use default settings with the `--local` alignment mode. We create indexes for four types of transcripts: miRNA, circRNA, piRNA, tRNA, and then for each of them:

1. We map the processed HiFi-Slide R2 reads `L2R2.trim_front_60.fastq`.
2. We select the uniquely mapped reads, i.e. reads "aligned exactly 1 time" in the Bowtie 2 log file. This is done by removing reads with auxiliary tag `XS`, i.e. reads that have other valid mappings.
3. We extract the fields of interest as HiFi-read identifier and transcript identifier.

```
for i in tRNA piRNA miRNA circRNA; do

bowtie2 \
-x bowtie2_index_$i \
-U L2R2.trim_front_60.fastq \
-S L2R2_$i"_mapped.sam" \
--un L2R2_$i"_unmapped.txt" \
--no-unal --threads 32 --local 2> L2R2_$i_mapped.log

samtools view L2R2_$i"_mapped.sam" | grep -v "XS:i:" > L2R2_$i"_uniquely_mapped.sam"

cut -f 1,3 L2R2_$i"_uniquely_mapped.sam" > L2R2_$i"_uniquely_mapped.txt"

done
```

We obtain four tab-separated txt files as output (as long as there are valid mapped HiFi-Slide reads for each of them). Example of output for piRNA:
```
MN00185:308:000H3YMVH:1:11101:11533:3933        piR-hsa-2277753
MN00185:308:000H3YMVH:1:11101:7721:4462         piR-hsa-2277753
MN00185:308:000H3YMVH:1:11101:20516:4919        piR-hsa-2277753
MN00185:308:000H3YMVH:1:11101:8487:6158         piR-hsa-2229595
```


## 3. Align HiFi-Slide R1 reads to deduplicated spatial barcodes

Aligner BWA is used to map HiFi-Slide R1 reads (L2R1) `L2R1.fastq` to deduplicated spatial barcodes `L1R1_dedup.fasta`. The minimum number of base match between L2R1 and L1R1 (spatial barcodes) `-k` is set at 80% of the length of spatial barcodes.

```
temp=$(less L1R1_dedup.fasta | head -2 | sed -n '2p')
min_base_match=$(echo "scale=4 ; ${#temp} * 0.8" | bc | awk '{printf("%.0f",$1)}')

bwa mem \
-a \
-k $min_base_match \
-t 32 \
L1R1_dedup L2R1.fastq > L2R1_L1R1_dedup.sam 2>L2R1_L1R1_dedup.log
```

In the case of a flowcell coming from MiniSeq, there is an additional step where we count the number of HiFi-Slide R1 reads mapped (flag 0 or 256) to spatial barcodes on the top (AAAL33WM5:1:1) and bottom surface (AAAL33WM5:1:2) and then select the surface with the highest number of aligned reads. This subsetted SAM file will go to the following steps.


## 4. Select HiFi-Slide R1 reads where R2 reads map over genes/transcripts

At this step, we use the HiFi-Slide read ID to select L2R1 reads aligned to spatial barcodes (L2R1_L1R1_dedup.sam) that have corresponding R2 ends (L2R2) mapped over genes/transcripts using the files generated at step 2. This step produces a filtered SAM file `L2R1_L1R1_dedup.filter.sam`.


## 5. Select HiFi-Slide R1 reads with highest alignment score
```
hifislida.pl L2R1_L1R1_dedup.filter.sam > L2R1_L1R1_dedup.hifislida.o 2>L2R1_L1R1_dedup.hifislida.e
```

**Arguments**  
1. Output SAM file from BWA.  

**Purpose**   
Read the SAM file output from BWA and collect the spatial barcodes aligned with each HiFi-Slide R1 reads at the highest alignment score. If a HiFi-Slide R1 read was aligned with N spatial barcodes that tied at the highest score, all the N spatial barcodes will be outputted except when N > 1,000. Spatial barcodes aligned with lower score would be discarded.

**Output**
Tab-separated file `L2R1_L1R1_dedup.hifislida.o` with the following columns:

- Column 1: HiFi-Slide read ID. 
- Column 2: Spatial barcode read ID. Each spatial barcode ID provide its spatial coordinates explicitly.  
- Column 3: Number of non-redundant spatial barcodes aligned at the highest score.  
- Column 4: Number of all the spatial barcodes aligned at the highest score.  
- Column 5: Highest alignment score between this HiFi-Slide R1 read and the aligned spatial barcode.


<!-- ## 4. Identify the Region of Interest (ROI)

### Rank the tiles by number of HiFi-Slide read pairs

```
hifislida2.pl L2R1_L1R1_dedup.hifislida.o L2R1_L1R1_dedup.sam > L2R1_L1R1_dedup.hifislida2.o 2>L2R1_L1R1_dedup.hifislida2.e
```

**Arguments**  
1. Output file produced by hifislida.pl.
2. Output file produced by BWA.

**Purpose**  
Count the number of univocally resolved HiFi reads per tile. A total of 6 X 11 tiles are available on NextSeq flowcell (sometimes they could be 6 X 14). We hypothesize that spatial barcodes on tiles covered by tissue should be mapped with more HiFi-Slide R1 reads than spatial barcodes outside the tissue covered region. To this end, we count the number of HiFi-Slide R1 reads per tile. To have a simplilified solution, we only consider HiFi-Slide R1 reads that have only one unique spatial barcode with the highest alignment score on the surface, which means considering only rows in `L2R1_L1R1_dedup.hifislida.o` where both columns 3 and 4 are equal to 1.       

**Output**  
Tab-separated file `L2R1_L1R1_dedup.hifislida2.o` with the following columns:

- Column 1: Identifier ">TILE".
- Column 2: Tile ID.
- Column 3: Number of spatially resolved HiFi-Slide read pairs per tile (ranked in descending order).
- Column 4: Number of mapped barcodes on the tile used for resolving the HiFi-Slide read pairs.

### Select tiles under ROI

For simplicity, we estimate the ROI as rectangular with a certain maximum and minimum size in terms of number of tiles per side. Next, we analyze all the possible ROI configurations between min and max size, and for each of them we perform a statistical test (Kolmogorov-Smirnov) between the number of HiFi-reads within and outside that ROI. ROIs with significant p-value are selected and sorted by avg_log2FC: the ROI with the highest avg_log2FC is selected and the tiles within it are given as output. If no significant ROI comes out, then all the tiles of the surface are considered as ROI.

```
select_tiles_in_ROI.r \
-i L2R1_L1R1_dedup.hifislida2.o \
-o ROI_tile_IDs.txt \
-f NextSeq \
--surface 1 \
--max_size_ROI 4 \
--min_size_ROI 2 \
--p_value 0.05
```

**Arguments**  
- `-i`: Input file, i.e. output from `hifislida2.pl`, tiles ranked by the number of spatially resolved HiFi-Slide read pairs per tile.
- `-o`: Output file, tile IDs under ROI.
- `-f`: Flowcell type, either MiniSeq, NextSeq.
- `--surface`: Flowcell surface where the tissue is. For NextSeq is `1` for MiniSeq it depends on which surface is the highest number of mapped HiFi-Slide R1 reads to barcodes. This parameter is automatically computed in the script and passed to this function.
- `--max_size_ROI`: the maximum size of a side of the ROI. For example, `max_size_ROI = 4` indicates that ROI can be at most 4 x 4 tiles.
- `--min_size_ROI`: the minimum size of a side of the ROI. For example, `min_size_ROI = 2` indicates that ROI shall be at least 2 x 2 tiles.
- `--p_value`: The p-value threshold to be used for statistical significance.
 -->

## 6. Match HiFi-Slide read pairs with spatial location 

```
hifislida3.pl \
L2R1_L1R1_dedup.hifislida.o \
L1R1_dup.txt > L2R1_L1R1.hifislida3.o
```

**Arguments**  
1. Output file produced by hifislida.pl.
2. The second output file from `surfdedup`. This file provides duplicate spatial barcodes whose IDs were not shown in the output from `hifislida.pl`.

**Purpose**  
With `hifislida.pl` we obtained aligned deduplicated spatial barcodes with the highest score for each HiFi-Slide read pair. These HiFi-Slide read pairs were considered as spatially resolved. Next, we use `hifislida3.pl` to print out the coordinates of all the aligned spatial barcodes on each tile. Notably, when multiple spatial barcodes share the same sequence but from different coordinates, `hifislida3.pl` prints out all their coordinates.

**Output**  
Tab-separated file `L2R1_L1R1.hifislida3.o` with the following columns:

- Column 1: HiFi-Slide read ID.
- Column 2: Tile ID (only tiles under the ROI provided as input).
- Column 3: X-coord on the tile (columns).
- Column 4: Y-coord on the tile (rows).
- Column 5: N as the number of total spatial coordinates where this HiFi-Slide read could be aligned to spatial barcodes. This is used to "weight" HiFi-Slide reads. For example, if a HiFi-Slide read has N = 8, it would be counted as 1/8 at any of these 8 coordinates.

Example:
```
HiFi_read_id    tile_id col     row     N
MN00185:308:000H3YMVH:1:21104:22534:6585        1308    20557   77594   4
MN00185:308:000H3YMVH:1:21104:22534:6585        1308    7797    22530   4
MN00185:308:000H3YMVH:1:21104:22534:6585        1209    30647   54019   4
MN00185:308:000H3YMVH:1:21104:22534:6585        1308    18891   36353   4
```


## 8. Integrate spatial coordinates and RNA information

As a final step, we integrate the outcomes from the sections above to associate a spatial coordinate with each expressed gene/transcript. This is done by performing a `join` of the output tables above using the HiFi-Slide read identifiers (after having sorted them by HiFi-Slide read ID). First, we generate tables where each HiFi-Slide read ID is associated with a spatial coordinate (coming from step 1) and a gene/transcript (coming from step 2). Finally, we aggregate those tables to obtain a list of unique "spot-gene" pairs with the corresponding gene expression levels.

The file with HiFi-Slide read spatial coordinates is `L2R1_L1R1.hifislida3.o`, which is sorted first to obtain `L2R1_L1R1.hifislida3.sort.o` which is used below.

```
cat L2R1_L1R1.hifislida3.o | sort -k 1 --parallel=32 -S 20G > L2R1_L1R1.hifislida3.sort.o
```

### Genome

```
cat HiFi_L2R2_genome.bed | sort -k 4 --parallel=32 -S 20G > HiFi_L2R2_genome.sort.bed

join -1 1 -2 1 -t $'\t' L2R1_L1R1.hifislida3.sort.o HiFi_L2R2_genome.sort.bed | cut -f 2,3,4,5,6,7,8 > HiFi_L2R2_genome_spatial.txt
```

Without the `cut` the first column would be the HiFi-Slide read ID (field used for `join`), however since we do not need this information anymore we prefer to remove that field before writing to output file to reduce memory occupation.

Tab-separated txt file with the following columns:

- Column 1: Tile ID.
- Column 2: X-coord on the tile (columns).
- Column 3: Y-coord on the tile (rows).
- Column 4: N as the number of total spatial coordinates where this HiFi-Slide read could be aligned to spatial barcodes.
- Column 5: Gene ID.
- Column 6: Gene name.
- Column 7: Gene biotype.

Example:

```
1109    45545   10979   2       ENSG00000230876.8       LINC00486       lncRNA
1109    45545   10979   4       ENSG00000230876.8       LINC00486       lncRNA
1109    13609   68183   1       ENSG00000230876.8       LINC00486       lncRNA
1109    14271   75095   22      ENSG00000230876.8       LINC00486       lncRNA
```

We then aggregate `HiFi_L2R2_genome_spatial.txt` to get unique (spot_i, gene_j) rows, and we calculate the expression level (counts) of gene_j in spot_i in this way: first we calculate `1/N` values in every row of `HiFi_L2R2_genome_spatial.txt`, then we sum the 1/N values corresponding to each (spot_i, gene_j) pair.

```
awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' HiFi_L2R2_genome_spatial.txt > HiFi_L2R2_genome_spatial.final.txt
```

Tab-separated txt file `HiFi_L2R2_genome_spatial.final.txt` with the following columns:

- Column 1: Tile ID.
- Column 2: X-coord on the tile (columns).
- Column 3: Y-coord on the tile (rows).
- Column 4: Gene ID.
- Column 5: Gene name.
- Column 6: Gene biotype.
- Column 7: Gene expression level (counts).

Example:

```
1109    45545   10979   ENSG00000230876.8       LINC00486       lncRNA  0.75
1109    13609   68183   ENSG00000230876.8       LINC00486       lncRNA  1
1109    14271   75095   ENSG00000230876.8       LINC00486       lncRNA  0.04545455
```

### Transcriptome

```
for i in tRNA piRNA miRNA circRNA; do

cat L2R2_$i"_uniquely_mapped.txt" | sort -k 1 --parallel=32 -S 20G > L2R2_$i"_uniquely_mapped.sort.txt"

join -1 1 -2 1 -t $'\t' L2R1_L1R1.hifislida3.sort.o L2R2_$i"_uniquely_mapped.sort.txt" > HiFi_L2R2_$i"_spatial.txt"

done
```

Tab-separated txt file with the following columns:

- Column 1: Tile ID (only tiles under the ROI provided as input).
- Column 2: X-coord on the tile (columns).
- Column 3: Y-coord on the tile (rows).
- Column 4: N as the number of total spatial coordinates where this HiFi-Slide read could be aligned to spatial barcodes.
- Column 5: Transcript identifier.

Example for piRNA:
```
1109	29549	34175	15	piR-hsa-2229595
1109	29549	34175	20	piR-hsa-2229595
1109	47135	43472	15	piR-hsa-2229595
1208	16960	53394	10	piR-hsa-2229595
```

Next, we perform the same aggregation method as above for transcripts.

```
for i in tRNA piRNA miRNA circRNA; do

awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$5]+=1/$4} END { for (i in array) {print i"\t" array[i]}}' HiFi_L2R2_$i"_spatial.txt" > HiFi_L2R2_$i"_spatial.final.txt"

done
```

Example for piRNA:
```
1109	29549	34175	piR-hsa-2229595 0.1166667
1109	47135	43472	piR-hsa-2229595 0.06666667
1208	16960	53394	piR-hsa-2229595 0.1
```



## QC metrics

All the QC metrics listed in section 4 [here](https://docs.google.com/document/d/1MvXPgTVzzeAEnmRXDRuaJMY-U_ENorQMPzqpLVOJWA0/edit?pli=1#) are saved in tab-separated txt files named respectively:

- `flowcell_ID.QC_metrics.txt`.
- `sample_name.QC_metrics.txt`.







