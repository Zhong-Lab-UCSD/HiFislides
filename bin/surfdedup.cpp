#include <fstream> 
#include <iostream> 
#include <sstream>
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
	unordered_map<string,string> read2sign;
	unordered_map<string,string> read2qscore;
	int nraw = 0;

	string surf;
	surf.assign(argv[1]);

	for (int i = 2; i < argc; i++) {
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

				stringstream sstream1;
				sstream1 << seq_id_2;
  		  		string seq_id_3;
    				std::getline(sstream1, seq_id_3, ' ');
				
				// seq_id_3.find(surf);
				std::size_t found = seq_id_3.find(surf);

				if(found!=std::string::npos) {
					
					nraw = nraw + 1;
					// readid2seq.insert( std::pair<string,string>(seq_id,raw_sequence) );
					auto autoit = read2freq.find(raw_sequence);
					if(autoit != read2freq.end()) {
						// exists
						int n = autoit->second;
						n = n + 1;
						autoit->second = n;

						auto it_read2id = read2id.find(raw_sequence);
						cerr << it_read2id->second << "\t" << seq_id_3 << "\n";

					} else {
						// not exists
						int i = 1;
						read2freq.insert( std::pair<string,int>(raw_sequence,i) );
						read2id.insert( std::pair<string,string>(raw_sequence,seq_id_3) );
					}
				}
			}
		gzclose(input_file);
		cerr << "Done " << fin << "\n";
	}
	int n1 = 0;
	unordered_map<string, int>::iterator it_read2freq = read2freq.begin();
	while(it_read2freq != read2freq.end()) {
		if(it_read2freq->second == 1) {
			n1 = n1 + 1;
		}
		auto it_read2id = read2id.find(it_read2freq->first);
		cout << ">" << it_read2id->second << "_" << it_read2freq->second << "\n";
		cout << it_read2freq->first;
		++it_read2freq;
	}
	cerr << "Raw reads: " << nraw << "\n";
	cerr << "Unique reads: " << n1 << "\n";
}
