# HiFislides
/*

Works by Stanley, Ekko, Pie, Brianne and many others in Zhong Lab

Hello world!

*/

Materials and Methods:

Aligner: STAR v2.7.5
two parameters use non-default setting for a less stringent criteria: --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2
(default for both are 0.66)

Genome and junction annotation: Ensembl release 104, Homo sapiens


# compile the c++ source code count_read_occurrence_in_Library_1.cpp:
# -lz is important for using zlib which allow c++ to read compressed fastq gz files
g++ -o count_read_occurrence_in_Library_1 count_read_occurrence_in_Library_1.cpp -lz
cp ./count_read_occurrence_in_Library_1 ~/bin
