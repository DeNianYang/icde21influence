#ifndef FUNCTION_H
#define FUNCTION_H

#include<cstdlib>
#include<iostream>
#include<string>
#include<fstream>
#include<queue>
#include<ctime>
#include<vector>
#include<map>
#include "Data.h"
using namespace std;

int str2int(string l);
double str2double(string l);
double get_random();

// read file
void read_all_file(map<int,Node*>& allNode,map<int,double>& allItem,double***& Sim,int mc,int ms,
					string ItemFile,string GraphFile,string metacFile,string metasFile,
					string PprefFile,string PextFile,string SimFile,string costFile);
void read_item(map<int,double>& allItem,string ItemFile);
void read_graph(map<int,Node*>& allNode,string GraphFile,int nitem,int mc,int ms);
void read_W_meta_c(map<int,Node*>& allNode,string metacFile);
void read_W_meta_s(map<int,Node*>& allNode,string metasFile);
void read_Ppref(map<int,Node*>& allNode,string PprefFile);
void read_Pext(map<int,Node*>& allNode,string PextFile);
void read_Seed(map<int,Node*>& allNode,vector<Seed*>& S,string SeedFile);
void read_Sim(double*** Sim,string SimFile);
void read_cost(map<int,Node*>& allNode,string costFile);

// diffusion
double diffusion(vector<Seed*>& S,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext,int mode);
void exertInf(Node* u,int x,vector<pair<Node*,int> >& newInf,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int mc,int ms,double eta,double Pminpref,double Pminext,int mode);
void UpdatePpref(Node* v,int x,double*** R,int nit);
void UpdateWeight(Node* v,int x,vector<pair<int,char> >& EA_pos,vector<pair<int,char> >& EA_neg,
					double*** Sim,int nit,int mc,int ms,double eta);

void rel(Node* v,double*** R,map<int,Node*>& allNode,double*** Sim,int nit,int mc,int ms);

void union_Inf(vector<pair<Node*,int> >& A,Node* v,int x);
void union_EA(vector<pair<int,char> >& A,int y,char c);
bool exist_item(vector<int>& A,int a);
void convert_seed(vector<Seed*>& Seed,vector<vector<pair<Node*,int> > >& S);
void InsertionSort(vector<Seed*>& A,int len);
void QuickSort(vector<Seed*>& A,int front,int end);
int Partition(vector<Seed*>& A,int front,int end);

// greedy
void Greedy(vector<Seed*>& GreedyS,map<int,Node*>& allNode,map<int,double>& allItem,double*** Sim,
			int N,int mc,int ms,double eta,double Pminpref,double Pminext,double b,int U_mode,fstream& f);
void CostGreedy(vector<Seed*>& S1,vector<Seed*>& U,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext,double b);
void USM(vector<Seed*>& S1_usm,vector<Seed*>& S1, map<int,Node*>& allNode,map<int,double>& allItem,
			double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext);
void SingleGreedy(vector<Seed*>& S3,vector<Seed*>& U,map<int,Node*>& allNode,map<int,double>& allItem,
					double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext,double b);


void max_greedy(vector<Seed*>& GreedyS,vector<Seed*>& S1,vector<Seed*>& S1_usm,vector<Seed*>& S2,vector<Seed*>& S3,
				map<int,Node*>& allNode,map<int,double>& allItem,double*** Sim,
				int N,int mc,int ms,double eta,double Pminpref,double Pminext);

// assignment
void assignment(vector<Seed*>& Smax,vector<Seed*>& S,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int nit,int mc,int ms,double eta,double Pminpref,double Pminext,
				int d,double theta,int T,fstream& f);
void TMI(vector<Seed*>& S,vector<vector<TM*> >& PTMs,map<int,Node*>& allNode,
			double*** Sim,int nit,int mc,int ms,int d,double theta,int T);
void UserClustering(vector<vector<Seed*> >& UC,vector<Seed*>& S,map<int,Node*>& allNode,int d);
void ItemClustering(vector<vector<Seed*> >& IC,vector<Seed*>& uc,double*** R);
int DRE(vector<Seed*>& S1,TM* tm,map<int,Node*>& allNode,map<int,double>& allItem,
		double*** Sim,int N,int nit,int mc,int ms,double eta,double Pminpref,double Pminext,int d);
double ImpactSample(int x,double*** R,char dir,map<int,double>& allItem,int N,int nit,int d);
void TDSI(vector<Seed*>& C,vector<Seed*>& S1,vector<Seed*>& last_S,TM* tm,map<int,Node*>& allNode,map<int,double>& allItem,
			double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext,int T,int min_T,int max_T);
pair<double,double> InfPrefSample(vector<Seed*>& Seeds,vector<Node*>& vtm,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext);
int distance(map<int,Node*>& allNode,Node* u,Node* v);
void distance(Node* u);
void RelAvg(double*** R,map<int,Node*>& allNode,double*** Sim,int nit,int mc,int ms);
double get_ratio(TM* a,TM* b);
void DFS(vector<TM*>& G,TM* n);
void SortTM(vector<TM*>& CC);
void union_Seed(vector<Seed*>& S,vector<Seed*>& C);
void union_Seed(vector<Seed*>& S,Seed* c);
void union_Node(vector<Node*>& S,vector<Node*>& C);
void union_Node(vector<Node*>& S,Node* c);
void union_vec_Seed(vector<vector<Seed*> >& UC,vector<Seed*>& C);
void diff_vec_Seed(vector<vector<Seed*> >& UC,vector<Seed*>& C);
void diff_Seed(vector<Seed*>& S,vector<Seed*>& C);
void diff_Seed(vector<Seed*>& S,Seed* c);
void union_Item(vector<int>& next,int z);
void UpdateRelAvg(vector<Seed*>& Seeds,vector<Node*>& VTM,double*** R,map<int,Node*>& allNode,map<int,double>& allItem,
					double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext);

// opt
void opt(vector<Seed*>& S,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int mc,int ms,double eta,double b,int T,fstream& f,int NUM_THREAD);
void EnumerateRounds(vector<vector<Seed*> >& next,vector<Seed*>& S1,Node* u,int x,int T);
bool seed_in_vec(Node* u,int x,int t,vector<Seed*>& S1);
void test_thread(void* p);

#endif
