
# this step read in raw reads from recycled flow cell and add "_N" to the identifier of each read.
# for a read with "_N", its sequence could be found N times on the surface.
# only reads whose identifier contain a string $surf would be used to determine N.
# by i = 1 and j = 1 we only raw reads from top surface of the lane 1.
#
surfdedup $surf *_L00$i\_R1_001.fastq.gz > L1R1Dedup_$i$j.fasta
