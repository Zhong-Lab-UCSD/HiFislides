#include <fstream> 
#include <iostream> 
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set> 
#include <utility>
#include "/home/linpei/bin/zlib.h"
#include <string>
#include <cmath> 
#include <stdio.h>
#include <iterator>
#include <unistd.h>
#include <sysexits.h>
#include <vector>
#include <algorithm>

# Reference: 
# https://doi.org/10.1093/bioinformatics/btaa112
# https://genome.cshlp.org/content/30/9/1364

# codes from https://github.com/daihang16/nubeamdedup
# were used and revised to print out the number of occurrence of each raw read in L1 library.


using namespace std; 

void operation(double prod[][2], bool yes) 
{
	if (yes) {
		prod[0][1] += prod[0][0]; 
		prod[1][1] += prod[1][0]; 
	} else {
		prod[0][0] += prod[0][1]; 
		prod[1][0] += prod[1][1]; 
	}
}

void matrix_multiplication_helper(double prod[][2], char * raw_seq, char nucleotide) {
	while (*raw_seq != '\n') {
		operation(prod, *raw_seq == nucleotide);
		raw_seq++;
	}
}

int main(int argc,char ** argv)
{
	unordered_set<double> vdat;
	unordered_map<double, string> nub2read;
	unordered_map<string, int> read2freq;
	vector<string> raw_read_id;
	cout << "Hello! Start Nubeamization!\n"; 
	for (int i = 1; i < argc; i++) {
		string fin;
		fin.assign(argv[i]);
		gzFile input_file = gzopen(fin.c_str(), "r"); 
		if (input_file == NULL) { printf("can't open %s file to read. \n", fin.c_str()); exit(EX_NOINPUT); }
		char * seq_id = new char[310];  
		char * raw_sequence = new char[310]; 
		char * seq_id_repeat = new char[310];
		char * tscore = new char[310];  
		double read_identifier = 0;
		unsigned int line_num = 0; // number of reads
		unsigned int n_uniq_reads = 0; // number of unique reads
		while (gzgets(input_file, seq_id, 310) &&
			gzgets(input_file, raw_sequence, 310) &&
			gzgets(input_file, seq_id_repeat, 310) &&
			gzgets(input_file, tscore, 310)) {
				double prod[2][2] = {{1,0},{0,1}};
				matrix_multiplication_helper(prod, raw_sequence, 'A');
				matrix_multiplication_helper(prod, raw_sequence, 'T');
				matrix_multiplication_helper(prod, raw_sequence, 'C');
				matrix_multiplication_helper(prod, raw_sequence, 'G');
				read_identifier = prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1];
				auto autoit = nub2read.find(read_identifier);
				if(autoit != nub2read.end()) {
					// exists
					int n = read2freq[autoit->second];
					n = n + 1;
					auto autoit2 = read2freq.find(autoit->second);
					autoit2->second = n;
				} else {
					// not exists
					nub2read.insert( std::pair<double,string>(read_identifier,seq_id) );
					int i = 1;
					read2freq.insert( std::pair<string,int>(seq_id,i) );
				}
		}
		gzclose(input_file);
	}
	unordered_map<string, int>::iterator it1 = read2freq.begin();
	while(it1 != read2freq.end()) {
		cout << it1->second << "\t" << it1->first;
		++it1;
	}
	cout << "Done! Happy Nubeamization!\n";	
}
