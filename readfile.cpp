#include<cstdlib>
#include<iostream>
#include<string>
#include<fstream>
#include<queue>
#include<ctime>
#include<vector>
#include<map>
#include "function.h"
using namespace std;

void read_all_file(map<int,Node*>& allNode,map<int,double>& allItem,double***& Sim,int mc,int ms,
					string ItemFile,string GraphFile,string metacFile,string metasFile,
					string PprefFile,string PextFile,string SimFile,string costFile)
{
	read_item(allItem,ItemFile); // read weight W_x
	read_graph(allNode,GraphFile,int(allItem.size()),mc,ms); // read graph and initialize parameters
	read_W_meta_c(allNode,metacFile);
	read_W_meta_s(allNode,metasFile);
	read_Ppref(allNode,PprefFile);
	read_Pext(allNode,PextFile);

	// create three dimension array Sim
	int nit = int(allItem.size());
	Sim = new double**[nit];
	for(int i=0;i<nit;i++)
	{
		Sim[i] = new double*[nit];
		for(int j=0;j<nit;j++)
			Sim[i][j] = new double[mc+ms];
	}
	read_Sim(Sim,SimFile);
	read_cost(allNode,costFile);

	// backup for reset
	map<int,Node*>::iterator it;
	for(it = allNode.begin();it != allNode.end();it++)
		it->second->backup(nit,mc,ms);
}

void read_item(map<int,double>& allItem,string ItemFile)
{
	fstream f;
	f.open(ItemFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));

	string s = "";
	int n = 0;
	double v = 0.0;
	while(!f.eof())
	{
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',')
			{
				n = str2int(s);
				s = "";
			}
			else
				s += buffer[i];
		}
		v = str2double(s);
		s = "";
		allItem[n] = v;
		f.getline(buffer,sizeof(buffer));
	}
}

void read_graph(map<int,Node*>& allNode,string GraphFile,int nitem,int mc,int ms)
{
	fstream f;
	f.open(GraphFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));

	string s = "";
	int u1 = 0;
	int u2 = 0;
	double v = 0.0;
	int count = 0;
	while(!f.eof())
	{
		count = 0;
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',')
			{
				if(count == 0)
					u1 = str2int(s);
				else if(count == 1)
					u2 = str2int(s);
				s = "";
				count++;
			}
			else
				s += buffer[i];
		}
		v = str2double(s);
		s = "";

		// insert edge
		Node* n1;
		Node* n2;
		if(allNode.find(u1) == allNode.end())
		{
			n1 = new Node(u1,nitem,mc,ms);
			allNode[u1] = n1;
		}
		else
			n1 = allNode[u1];
		if(allNode.find(u2) == allNode.end())
		{
			n2 = new Node(u2,nitem,mc,ms);
			allNode[u2] = n2;
		}
		else
			n2 = allNode[u2];
		n1->nout[n2] = v;
		n2->nin.push_back(n1);
		// insert edge end

		f.getline(buffer,sizeof(buffer));
	}
}

void read_W_meta_c(map<int,Node*>& allNode,string metacFile)
{
	fstream f;
	f.open(metacFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));

	string s = "";
	int n = 0;
	int m = 0;
	double v = 0.0;
	int count = 0;
	while(!f.eof())
	{
		count = 0;
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',')
			{
				if(count == 0)
					n = str2int(s);
				else if(count == 1)
					m = str2int(s);
				s = "";
				count++;
			}
			else
				s += buffer[i];
		}
		v = str2double(s);
		s = "";
		allNode[n]->W_meta_c[m] = v;

		f.getline(buffer,sizeof(buffer));
	}
}

void read_W_meta_s(map<int,Node*>& allNode,string metasFile)
{
	fstream f;
	f.open(metasFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));

	string s = "";
	int n = 0;
	int m = 0;
	double v = 0.0;
	int count = 0;
	while(!f.eof())
	{
		count = 0;
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',')
			{
				if(count == 0)
					n = str2int(s);
				else if(count == 1)
					m = str2int(s);
				s = "";
				count++;
			}
			else
				s += buffer[i];
		}
		v = str2double(s);
		s = "";
		allNode[n]->W_meta_s[m] = v;

		f.getline(buffer,sizeof(buffer));
	}
}

void read_Ppref(map<int,Node*>& allNode,string PprefFile)
{
	fstream f;
	f.open(PprefFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));

	string s = "";
	int user = 0;
	int item = 0;
	double v = 0.0;
	int count = 0;
	while(!f.eof())
	{
		count = 0;
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',')
			{
				if(count == 0)
					user = str2int(s);
				else if(count == 1)
					item = str2int(s);
				s = "";
				count++;
			}
			else
				s += buffer[i];
		}
		v = str2double(s);
		s = "";
		allNode[user]->Ppref[item] = v;

		f.getline(buffer,sizeof(buffer));
	}
}

void read_Pext(map<int,Node*>& allNode,string PextFile)
{
	fstream f;
	f.open(PextFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));

	string s = "";
	int user = 0;
	int item1 = 0;
	int item2 = 0;
	double v = 0.0;
	int count = 0;
	while(!f.eof())
	{
		count = 0;
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',')
			{
				if(count == 0)
					user = str2int(s);
				else if(count == 1)
					item1 = str2int(s);
				else if(count == 2)
					item2 = str2int(s);
				s = "";
				count++;
			}
			else
				s += buffer[i];
		}
		v = str2double(s);
		s = "";
		allNode[user]->Pext[item1][item2] = v;

		f.getline(buffer,sizeof(buffer));
	}
}

void read_Seed(map<int,Node*>& allNode,vector<Seed*>& S,string SeedFile)
{
	fstream f;
	f.open(SeedFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));

	string s = "";
	int user = 0;
	int item = 0;
	int seed = 0.0;
	int count = 0;
	//int seednum = -1;
	while(!f.eof())
	{
		count = 0;
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',')
			{
				if(count == 0)
					user = str2int(s);
				else if(count == 1)
					item = str2int(s);
				s = "";
				count++;
			}
			else
				s += buffer[i];
		}
		seed = str2int(s);
		s = "";
		/*if(seed > seednum)
		{
			vector<pair<Node*,int> > x;
			S.push_back(x);
			seednum = seed;
		}*/
		Seed* newseed = new Seed(allNode[user],item,seed);
		S.push_back(newseed);
		//pair<Node*,int> p(allNode[user],item);
		//S[seed].push_back(p);

		f.getline(buffer,sizeof(buffer));
	}
}

void read_Sim(double*** Sim,string SimFile)
{
	fstream f;
	f.open(SimFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));

	string s = "";
	int item1 = 0;
	int item2 = 0;
	int meta = 0;
	double v = 0.0;
	int count = 0;
	while(!f.eof())
	{
		count = 0;
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',')
			{
				if(count == 0)
					item1 = str2int(s);
				else if(count == 1)
					item2 = str2int(s);
				else if(count == 2)
					meta = str2int(s);
				s = "";
				count++;
			}
			else
				s += buffer[i];
		}
		v = str2double(s);
		s = "";
		Sim[item1][item2][meta] = v;

		f.getline(buffer,sizeof(buffer));
	}
}

void read_cost(map<int,Node*>& allNode,string costFile)
{
	fstream f;
	f.open(costFile.c_str(),ios::in);
	char buffer[1000];
	f.getline(buffer,sizeof(buffer));

	string s = "";
	int user = 0;
	int item = 0;
	double v = 0.0;
	int count = 0;
	while(!f.eof())
	{
		count = 0;
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',')
			{
				if(count == 0)
					user = str2int(s);
				else if(count == 1)
					item = str2int(s);
				s = "";
				count++;
			}
			else
				s += buffer[i];
		}
		v = str2double(s);
		s = "";
		allNode[user]->cost[item] = v;

		f.getline(buffer,sizeof(buffer));
	}
}
