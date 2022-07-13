#include<cstdlib>
#include<iostream>
#include<string>
#include<fstream>
#include<queue>
#include<ctime>
#include<vector>
#include<map>
#include "Data.h"
#include "function.h"
using namespace std;

int main(int argc,char* argv[])
{
	clock_t time1 = clock();
	// random seed
	srand(time(NULL));
	string paramFile = argv[1];

	// Basic info of user and item
	string ItemFile; // include weight
	string GraphFile; // include fanin fanout
	string metacFile;
	string metasFile;
	string PprefFile;
	string PextFile;
	string SimFile;
	int N = -1;
	double eta = -1.0;
	int mc = -1;
	int ms = -1;
	int U_mode = -1;
	double Pminpref = -1.0;
	double Pminext = -1.0;
	double b = -1.0;
	int d = -1;
	double theta = -1;
	int T = -1;
	string costFile;
	string outputFile;

	fstream f;
	string s = "";
	f.open(paramFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));
	int line = 0;
	while(!f.eof())
	{
		line++;
		for(int i=0;buffer[i]!='\0';i++)
			s += buffer[i];
		switch(line)
		{
			case 1: ItemFile = s; break;
			case 2: GraphFile = s; break;
			case 3: metacFile = s; break;
			case 4: metasFile = s; break;
			case 5: PprefFile = s; break;
			case 6: PextFile = s; break;
			case 7: SimFile = s; break;
			case 8: costFile = s; break;
			case 9: outputFile = s; break;
			case 10: N = str2int(s); break;
			case 11: eta = str2double(s); break;
			case 12: mc = str2int(s); break;
			case 13: ms = str2int(s); break;
			case 14: Pminpref = str2double(s); break;
			case 15: Pminext = str2double(s); break;
			case 16: U_mode = str2int(s); break;
			case 17: b = str2double(s); break;
			case 18: d = str2int(s); break;
			case 19: theta = str2double(s); break;
			case 20: T = str2int(s); break;
		}
		s = "";
		f.getline(buffer,sizeof(buffer));
	}
	cout << "====================Parameters====================" << endl;
	cout << "N = " << N << ", eta = " << eta << ", meta_c_n = " << mc << ", meta_s_n = " << ms << endl;
	cout << "Pminpref = " << Pminpref << ", Pminext = " << Pminext << ", b = " << b << ", d = " << d << endl;
	cout << "theta = " << theta << ", T = " << T << endl;
	cout << "==================================================" << endl;
	if(U_mode == 0)
		cout << "Compute S1, S3" << endl;
	else
		cout << "Compute S1, S1', S2, S3" << endl;
	f.close();

	// Reading Files
	cout << "Reading Files ..." << endl;
	map<int,Node*> allNode;
	map<int,double> allItem; // W_x
	double*** Sim;
	read_all_file(allNode,allItem,Sim,mc,ms,ItemFile,GraphFile,metacFile,metasFile,PprefFile,PextFile,SimFile,costFile);

	int nit = int(allItem.size());
	fstream outputfile;
	outputfile.open(outputFile.c_str(),ios::out);

	// Greedy
	cout << "Greedy Algorithm ..." << endl;
	vector<Seed*> GreedyS;
	Greedy(GreedyS,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,b,U_mode,outputfile);

	// assignment
	cout << "Assignment ..." << endl;
	vector<Seed*> Smax;
	assignment(Smax,GreedyS,allNode,allItem,Sim,N,nit,mc,ms,eta,Pminpref,Pminext,d,theta,T,outputfile);

	outputfile.close();

	cout << "Time Usage: " << double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds." << endl;
}

