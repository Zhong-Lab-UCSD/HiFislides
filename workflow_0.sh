#
# compile the c++ source code count_read_occurrence_in_Library_1.cpp:
#
# -lz is important for using zlib which allow c++ to read compressed fastq gz files
#
#
# g++ -o cro count_read_occurrence_in_Library_1.cpp -lz
# cp ./cro ~/bin

cro *_*_L00*_R1_001.fastq.gz > Library1_readundancy_check.o 2>Library1_readundancy_check.e

# then we used R to plot the distribution
