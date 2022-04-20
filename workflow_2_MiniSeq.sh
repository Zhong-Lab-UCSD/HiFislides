
L1R1=./lib1/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz 
L2R1=./lib2/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz

# gunzip -c $L1R1 | grep "MN00185" | cut -d " " -f 1 | cut -d ":" -f 5 | sort | uniq

nohup bwa index -p L1R1 $L1R1 > bwaindexo 2>bwaindexe &

k=54

date
bwa mem L1R1 $L2R1 -a -k $k -t 64 > L2R1toL1_k$k.sam 2>L2R1toL1_k$k.e;
date

grep -P "^MN00185:260:000H3W2LL:\d:\d+:\d+:\d+\t0\t" L2R1toL1_k$k.sam | cut -f 1,2,3 > L2R1toL1_k$k.hits
grep -P "^MN00185:260:000H3W2LL:\d:\d+:\d+:\d+\t256\t" L2R1toL1_k$k.sam | cut -f 1,2,3 >> L2R1toL1_k$k.hits

gethificoord.py L2R1toL1_k$k.hits > one_hit_hifi_reads_nspot2tiles.L 2>n_hits_hifi_reads.L
# sort -k 2nr one_hit_hifi_reads_nspot2tiles.L

################################################################################################################
date;
countreadundancy $L1R1 > Library1_readundancy_check.o 2>Library1_readundancy_check.e;
date
if [ -e unique_coord_per_tile.o ]; 
then 
echo "Found";
rm unique_coord_per_tile.o; 
fi

for j in `gunzip -c $L1R1 | grep "MN00185" | cut -d " " -f 1 | cut -d ":" -f 5 | sort | uniq`;
do
tile=FVV:1:$j:
n=`grep $tile Library1_readundancy_check.o | grep -P "\s1$" | wc -l` 
echo $tile $n >> unique_coord_per_tile.o
done
