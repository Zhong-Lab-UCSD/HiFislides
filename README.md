# HiFi-Slides Sequencing


# Works by Stanley, Ekko, Pie, Riccardo, Brianne and many others in Zhong Lab.


# variables
flowcell=AAAL33WM5  
i=1  
j=1  
surf=$flowcell:$i:$j   

# custom external scripts


**surfdedup**  
(1) input: 
string identifier of the surface  
fastqs of raw reads  
(2) Goal: this step read in raw reads from recycled flow cell and add "_N" to the identifier of each read. For a read whose sequence could be found N times on the surface, its identifier would be ended with "_N". Only reads whose identifier contain a string $surf would be used to determine N. By i = 1 and j = 1 we only raw reads from top surface of the lane 1. 


**finduniqread.pl**  
(1) input: the output fasta from surfdedup  
(2) Goal: this step will extract all raw reads with a "_1" in the original read identifier. 

**getgenefromgtf.pl**  

**hifia_asort.pl**  
(1) input:  output sam file from BWA and ID of the flowcell  
(2) Goal: read the sam file and output spatial barcode/coordinates that aligned to each HiFi R1 read.  
(3) Ouput:   
column 1 - spatial barcodes,  
column 2 - a useless number as placeholder,  
column 3 - identifiers of HiFi read pairs,  
column 4 - the number of aligned spatial coordinates for this HiFi R1 read  


**hifia_1n_marker_per_spot.pl**  
(1) input: 
(2) Goal: for each hifi read pair,hifia_1n_marker_per_spot.pl read the spatial barcodes aligned by R1, and genes mapped by R2. Then this script iterated on each spatial coordiante and print all related hifi reads and related genes, creating a spot-to-gene table.  
(3) Output:  
column 1 - ID of tile  
column 2 - X axis coordinate (columns)  
column 3 - Y axis coordinate (rows)  
column 4 - ID of spatial barcode  
column 5 - Ensembl Gene ID  
column 6 - official gene symbol  
column 7 - sum of read counts that mapped to this gene at this spot  
column 8 - ID of Hifi read pair. (when multiple HiFi read pairs could be mapped to 1 gene at the same coordiante, one of them would be shown randomly)    