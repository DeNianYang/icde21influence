#ifndef DATA_H
#define DATA_H

#include<cstdlib>
#include<iostream>
#include<string>
#include<fstream>
#include<queue>
#include<ctime>
#include<map>
#include<vector>
using namespace std;

class Node{
	public:
		Node(){};
		// ni: # of items
		// mc: # of meta_c
		// ms: # of meta_s
		// nn: # of nodes
		Node(int n,int ni,int mc,int ms)
		{
			num = n;
			W_meta_c = new double[mc]();
			W_meta_s = new double[ms]();
			Ppref = new double[ni]();
			Pext = new double*[ni];
			for(int i=0;i<ni;i++)
				Pext[i] = new double[ni]();
			active = new bool[ni]();
			for(int i=0;i<ni;i++)
				active[i] = false;
			cost = new double[ni]();

			B_W_meta_c = new double[mc]();
			B_W_meta_s = new double[ms]();
			B_Ppref = new double[ni]();
			dist = -1;
		}

		void backup(int ni,int mc,int ms)
		{
			for(int i=0;i<mc;i++)
				B_W_meta_c[i] = W_meta_c[i];
			for(int i=0;i<ms;i++)
				B_W_meta_s[i] = W_meta_s[i];
			for(int i=0;i<ni;i++)
				B_Ppref[i] = Ppref[i];
		}
		void reset(int ni,int mc,int ms,int mode)
		{
			if(mode == 0)
			{
				for(int i=0;i<mc;i++)
					W_meta_c[i] = B_W_meta_c[i];
				for(int i=0;i<ms;i++)
					W_meta_s[i] = B_W_meta_s[i];
				for(int i=0;i<ni;i++)
					Ppref[i] = B_Ppref[i];
			}
			for(int i=0;i<ni;i++)
				active[i] = false;
		}
	
		int num;
		double* 	W_meta_c;
		double* 	W_meta_s;
		double* 	Ppref;
		double** 	Pext;
		bool* 		active;
		double*	 	cost;

		double* 	B_W_meta_c;
		double* 	B_W_meta_s;
		double* 	B_Ppref;

		map<Node*,double> 	nout; // Pact
		vector<Node*> 		nin;
		int 				dist;
};

class Seed
{
	public:
		Seed(){}
		Seed(Node* n,int y,int a){u = n;x = y;t = a;}
		
		Node* u;
		int x;
		int t;
};

class TM
{
	public:
		TM(){}
		TM(vector<Seed*>& ic,int d)
		{
			visit = false;
			N_TM = ic;
			T_TM = 0;
			r_TM = 0.0;
			for(size_t i=0;i<ic.size();i++) // VTM <- {u|(u,x,t) in ic}
				V_TM.push_back(ic[i]->u);
			int h = 0;
			vector<Node*> now = V_TM; // now <- V_TM
			vector<Node*> next;
			map<Node*,double>::iterator it;
			while(h <= d && now.size() != 0) // while h <= d and now != {}
			{
				next.clear(); // next <- {}
				for(size_t i=0;i<now.size();i++) // for u in now do
				{
					Node* u = now[i];
					for(it = u->nout.begin();it != u->nout.end();it++) // for u's out neighbor v
					{
						// v is not in VTM U next
						Node* v = it->first;
						bool exist = false;
						for(size_t j=0;j<V_TM.size();j++)
						{
							if(v == V_TM[j])
							{exist = true;break;}
						}
						if(exist){continue;}
						for(size_t j=0;j<next.size();j++)
						{
							if(v == next[j])
							{exist = true;break;}
						}
						if(exist){continue;}

						// next <- next U {v}
						next.push_back(v);
					}
				}
				// VTM <- VTM U next
				for(size_t i=0;i<next.size();i++)
				{
					bool exist = false;
					for(size_t j=0;j<V_TM.size();j++)
					{
						if(next[i] == V_TM[j])
						{
							exist = true;
							break;
						}
					}
					if(!exist)
						V_TM.push_back(next[i]);
				}
				now = next;
				h++;
			}
		}
		
		vector<Seed*> 	N_TM;
		int 			T_TM;
		double 			r_TM;
		vector<Node*> 	V_TM;
		vector<TM*>		link;
		bool 			visit;
};
class Param
{
	public:
		Param(){};

		vector<Seed*>* U;
		map<int,Node*>* allNode;
		map<int,double>* allItem;
		vector<Seed*>* S1;
		vector<Seed*>* S;
		vector<vector<Seed*> >* now;
		double*** Sim;
		int N;
		int mc;
		int ms;
		double eta;
		int T;
		double b;
		double* max;
		int proc_num;
};

#endif
