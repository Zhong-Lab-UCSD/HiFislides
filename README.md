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

**hifia_1n_marker_per_spot.pl**  

