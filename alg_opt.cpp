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
	//double Pminpref = -1.0;
	//double Pminext = -1.0;
	double b = -1.0;
	//int d = -1;
	//double theta = -1;
	int T = -1;
	int NUM_THREAD = 1;
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
			case 14: b = str2double(s); break;
			case 15: T = str2int(s); break;
			case 16: NUM_THREAD = str2int(s); break;
		}
		s = "";
		f.getline(buffer,sizeof(buffer));
	}
	cout << "====================Parameters====================" << endl;
	cout << "N = " << N << ", eta = " << eta << ", meta_c_n = " << mc << ", meta_s_n = " << ms << endl;
	cout << ", b = " << b << ", T = " << T << endl;
	cout << "==================================================" << endl;
	f.close();

	// Reading Files
	cout << "Reading Files ..." << endl;
	map<int,Node*> allNode;
	map<int,double> allItem; // W_x
	double*** Sim;
	read_all_file(allNode,allItem,Sim,mc,ms,ItemFile,GraphFile,metacFile,metasFile,PprefFile,PextFile,SimFile,costFile);

	fstream outputfile;
	outputfile.open(outputFile.c_str(),ios::out);

	vector<Seed*> S;
	S.clear();
	opt(S,allNode,allItem,Sim,N,mc,ms,eta,b,T,outputfile,NUM_THREAD);
	//for(size_t i=0;i<S.size();i++)
	//	cout << S[i]->u->num << " " << S[i]->x << " " << S[i]->t << endl;
	outputfile.close();

	cout << "Time Usage: " << double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds." << endl;
}

