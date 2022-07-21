flowcell=AAAL33WM5
i=1
j=1
surf=$flowcell:$i:$j

cd /mnt/extraids/OceanStor-0/linpei/hifi/data_12/lib1/raw2
surfdedup $surf *_L00$i\_R1_001.fastq.gz > L1R1Dedup_$i$j.fasta
finduniqread.pl L1R1Dedup_$i$j.fasta > /mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/L1R1Uniq_$i$j.fasta

# R1 reads of HiFi sequencing - spatial end/coordiante end
L2R1=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
# R2 reads of HiFi sequencing - RNA end
L2R2=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/raw/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz

bwa index -p L2R1 $L2R1 > bwaindexL2R1o 2>bwaindexL2R1e

spatialbarcode=/mnt/extraids/OceanStor-0/linpei/hifi/data_14/lib2/L1R1Uniq_$i$j.fasta

seq=L1R1Uniq_11
k=40
sam=$seq\_L2R1_ak$k.sam
bwa mem L2R1 $spatialbarcode -a -k $k -t 64 > $sam 2>$seq\_L2R1_ak$k.same

mwd=/mnt/extraids/OceanStor-0/linpei
hg38=$mwd/imc/HG38
prefix=L2R2_000_

STAR --runThreadN 32 --genomeDir $hg38 --readFilesIn $L2R2 --quantMode GeneCounts --readFilesCommand zcat \
--outFileNamePrefix $prefix \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 > starlogo 2>starlogoe

grep "SN:" $prefix\Aligned.out.sam > $prefix\Aligned.NH1.sam
grep ":STAR" $prefix\Aligned.out.sam >> $prefix\Aligned.NH1.sam 
grep -P "NH:i:1\t" $prefix\Aligned.out.sam >> $prefix\Aligned.NH1.sam 

filetag=$prefix\Aligned.NH1
sam=$filetag.sam
bam=$filetag.bam
genicreadfile=$filetag\gene.L
samtools view -S -b $sam --threads 16 > $bam 2>>anye
gtf=$mwd/genome/release104/Homo_sapiens.GRCh38.104.gtf

getgenefromgtf.pl $gtf ENSG > genensg104.b 2>>anye

cat genensg104.b | perl -p -e "s/:/\t/g" | cut -f 1,2,3,4 > genensg104clean.b
bedtools intersect -a $bam -b genensg104clean.b -wb -bed > $genicreadfile
cut -f 1,2,3,4,16 $genicreadfile > $filetag\gene5columns.L

hifi2gene=L2R2_000_Aligned.NH1gene5columns.L

hifia_asort.pl $sam $flowcell > $seq\_L2R1_ak$k\_mappedspot_1n.L
ensgname=/mnt/extraids/OceanStor-0/linpei/genome/release104/ensg2name38104.txt
hifia_1n_marker_per_spot.pl $seq\_L2R1_ak$k\_mappedspot_1n.L $hifi2gene $flowcell hifi2gene.G $ensgname > Output_spot_to_gene.A

