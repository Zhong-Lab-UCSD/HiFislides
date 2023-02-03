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

using namespace std; 


int main(int argc,char ** argv)
{
	unordered_map<string, int> read2freq;
	unordered_map<string,string> read2id;
	for (int i = 1; i < argc; i++) {
		string fin;
		fin.assign(argv[i]);
		gzFile input_file = gzopen(fin.c_str(), "r"); 
		if (input_file == NULL) { printf("can't open %s file to read. \n", fin.c_str()); exit(EX_NOINPUT); }
		char * seq_id = new char[300];  
		char * raw_sequence = new char[300]; 
		char * seq_id_repeat = new char[300];
		char * tscore = new char[300];  
		double read_identifier = 0;
		unsigned int line_num = 0; // number of reads
		unsigned int n_uniq_reads = 0; // number of unique reads
		while (gzgets(input_file, seq_id, 300) &&
			gzgets(input_file, raw_sequence, 300) &&
			gzgets(input_file, seq_id_repeat, 300) &&
			gzgets(input_file, tscore, 300)) {
				std::string seq_id_2 = std::string(seq_id);
				seq_id_2.erase(seq_id_2.length()-1);
				seq_id_2.erase(0,1);
				// readid2seq.insert( std::pair<string,string>(seq_id,raw_sequence) );
				auto autoit = read2freq.find(raw_sequence);
				if(autoit != read2freq.end()) {
					// exists
					int n = autoit->second;
					n = n + 1;
					autoit->second = n;
				} else {
					// not exists
					int i = 1;
					read2freq.insert( std::pair<string,int>(raw_sequence,i) );
					read2id.insert( std::pair<string,string>(raw_sequence,seq_id_2) );
				}
		}
		gzclose(input_file);
	}
	unordered_map<string, int>::iterator it1 = read2freq.begin();
	while(it1 != read2freq.end()) {
		auto autoitt = read2id.find(it1->first);
		cout << autoitt->second << " " << it1->second << "\n";
		++it1;
	}
}
