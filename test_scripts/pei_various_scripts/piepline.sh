refdir=/mnt/extraids/OceanStor-0/linpei/genome/release104/STARhg38
r2=/mnt/extraids/OceanStor-0/linpei/hifi/data_1/15oct21_lib_2/L2R2.fastq 
op=L2R2test1
lib1sam=$mwd/hifi/data_1/15oct21_lib_1/library2_read1Aligned.out.sam

STAR --runThreadN 32 --genomeDir $refdir --readFilesIn $r2 --outFileNamePrefix $op --quantMode GeneCounts --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2
# grep -P "NH:i:1\s" L2R2test1Aligned.out.sam > L2R2test1Aligned.out.NHi1.sam

starout=L2R2test1ReadsPerGene.out.tab

grep -P "\tgene\t" $mwd/genome/release104/Homo_sapiens.GRCh38.104.gtf > $mwd/genome/release104/Homo_sapiens.GRCh38.104.gene.gtf
gtf=$mwd/genome/release104/Homo_sapiens.GRCh38.104.gene.gtf

# sam=L2R2test1Aligned.out.NHi1.sam
# determine which gene can be detected by Library 2 Read 2:
# nohup getgenefromsamfi.pl $starout $gtf $sam $j > stargene2library2read.$j.out 2>stargene2library2readcount.$j.out &
# getgenefromstarout.pl $starout $gtf

getgenefromstarout.pl $starout $gtf > locusdetectedbystar.o

# cut -f 1 stargene2library2.out | awk '!seen[$1]++'
# hifiplotjobg.pl $mwd/hifi/data_1/15oct21_lib_1/library2_read1Aligned.out.sam stargene2library2read.out > L2R1forgeneonmap.o
# hifiplotjobg.pl $mwd/hifi/data_1/15oct21_lib_1/library2_read1Aligned.out.sam stargene2library2read.out > L2R1forgeneonmap.o
# hifiplotjobg.pl $mwd/hifi/data_1/15oct21_lib_1/library2_read1Aligned.out.sam stargene2library2read.out o1 o2 2>e1

sam0=L2R2test1Aligned.out.sam
bam=L2R2test1Aligned.out.bam
sor=L2R2test1Aligned.out.sort.bam

samtools view -S -b $sam0 --threads 16 > $bam
samtools sort $bam -o $sor
samtools index $sor

samview.pl locusdetectedbystar.o $sor > locusdetectedbystar.sam 2>locusdetectedbystar.fullline.sam
hifiplotjobg.pl $mwd/hifi/data_1/15oct21_lib_1/library2_read1Aligned.out.sam locusdetectedbystar.sam o1x o2x > summary_11262021 2>e1x
