#include<iostream>
#include<fstream>
#include <math.h>
using namespace std;

int main(int argc,char ** argv)
{   

int row;
int col = 2;

row = atoi(argv[2]);
int** mat = new int*[row];

for(int i =0;i < row;i++) {
	mat[i] = new int[2];
}

double* dis;
dis  = new double[row];
for(int i = 0;i < row;i++) {
	dis[i] = 10000;
}

//Opening the file
string fin;
fin.assign(argv[1]);
ifstream inputfile(fin);


int start = atoi(argv[3]);
  int end = atoi(argv[4]);

// ifstream inputfile("2dinput.txt");
if (!inputfile.is_open()) 
	cout<<"Error opening file" ;
	//Defining the loop for getting input from the file
	for (int r = 0; r < row; r++) //Outer loop for rows
	{
		for (int c = 0; c < col; c++) //inner loop for columns
		{
			inputfile >> mat[r][c];  //Take input from file and put into mat
		}
	}
	for (int i = start; i < end; i++)
	{
		for (int j = 0; j < row; j++)
		{
			if(j != i) {
				double dx = abs(mat[i][0] - mat[j][0]);
				double dy = abs(mat[i][1] - mat[j][1]);
				if(dx <= 100 && dy <= 100) {
					double d = pow(dx,2) + pow(dy,2);
					double d1 = sqrt(d);
					if(d1 < dis[i]) {
						dis[i] = d1;
					}
					// cout << std::fixed << i << " " << j << " " << d << " " << d1 << "\n";
				}
			}
		}
	}
	for(int i = start;i < end;i++) {
		cout << i << " " << dis[i] << "\n";
	}
}
