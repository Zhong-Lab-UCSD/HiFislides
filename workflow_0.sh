#
# compile the c++ source code count_read_occurrence_in_Library_1.cpp:
#
# -lz is important for using zlib which allow c++ to read compressed fastq gz files
#
# g++ -o countreadundancy countreadundancy.cpp -lz
# cp ./countreadundancy ~/bin
# g++ -o countreadundancx countreadundancy_1.cpp -lz
# cp ./countreadundancx ~/bin

# L: lane
# T: tile

L=$1
T=$2

if [ $T eq 0 ]
then
	date;
	countreadundancy *_*_L00*_R1_001.fastq.gz > Library1_readundancy_check.o 2>Library1_readundancy_check.e
	date
else
	date
	tile=:$L:$T:
	countreadundancx $tile *_*_L00$L\_R1_001.fastq.gz > Library1_readundancy_checkLane$L\_Tile$T.o 2>Library1_readundancy_checkLane$L\_Tile$T.e
	date
fi

