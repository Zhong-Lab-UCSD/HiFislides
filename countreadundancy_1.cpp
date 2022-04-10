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
	std::string tile;
	tile = argv[1];

	int xmax = 0;
	int ymax = 0;

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
				stringstream sstream2;

				sstream1 << seq_id_2;
    				string seq_id_3;
    				std::getline(sstream1, seq_id_3, ' ');

				sstream2 << seq_id_3;
				std::string cluster[7] = { "machine","run","flowcell","Lane","tile","x","y" };

				for(int i=0;i < 7;i++) {
					std::getline(sstream2, cluster[i], ':');
				}
				int x_ = std::stoi(cluster[5]);
				int y_ = std::stoi(cluster[6]);
				if(x_ > xmax) {
					xmax = x_;
				}
				if(y_ > ymax) {
					ymax = y_;
				}
				// too slowly
				// std::regex self_regex(tile,std::regex_constants::ECMAScript | std::regex_constants::icase);
				// if (std::regex_search(seq_id_2,self_regex)) {
				std::size_t found = seq_id_2.find(tile);
  				if(found!=std::string::npos) {
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
					/*
					stringstream sstream1;
					stringstream sstream2;

					sstream1 << seq_id_2;
    					string seq_id_3;
    					std::getline(sstream1, seq_id_3, ' ');
					
					sstream2 << seq_id_3;

					std::string cluster[7] = { "machine","run","flowcell","Lane","tile","x","y" };
					for(int i=0;i < 7;i++) {
						std::getline(sstream2, cluster[i], ':');
					}
					
					for(int i=0;i < 7;i++) {
						cout << cluster[i] << ",";
					}
					cout << "\n";
					*/
				}
			}
		gzclose(input_file);
	}

	cerr << "Tile " << tile << " " << xmax << " " << ymax << "\n";

	unordered_map<string, int>::iterator it1 = read2freq.begin();
	while(it1 != read2freq.end()) {
		auto autoitt = read2id.find(it1->first);
		string seq_id = autoitt->second;

		stringstream sstream1;
		stringstream sstream2;

		sstream1 << seq_id;
    		string seq_id_3;
    		std::getline(sstream1, seq_id_3, ' ');
					
		sstream2 << seq_id_3;

		std::string cluster[7] = { "machine","run","flowcell","Lane","tile","x","y" };
		for(int i=0;i < 7;i++) {
			std::getline(sstream2, cluster[i], ':');
		}
		cout << seq_id_3 << " " << cluster[4] << " " << cluster[5] << " " << cluster[6] << " " << "Times:" << " " << it1->second << "\n";
		++it1;
	}
}
