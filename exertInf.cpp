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

void exertInf(Node* u,int x,vector<pair<Node*,int> >& newInf,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int mc,int ms,double eta,double Pminpref,double Pminext,int mode)
{
	// newInf <- {}
	newInf.clear();
	vector<pair<int,char> > EA_pos;
	vector<pair<int,char> > EA_neg;

	map<Node*,double>::iterator it;
	map<int,double>::iterator it2;
	int nit = int(allItem.size());

	double*** R = new double**[nit];
	for(int i=0;i<nit;i++)
	{
		R[i] = new double*[nit];
		for(int j=0;j<nit;j++)
			R[i][j] = new double[2];
	}

	double rd = 0.0;
	for(it = u->nout.begin();it != u->nout.end();it++) // for v in N_out(u) do
	{
		Node* v = it->first;

	    rd = double(rand())/double(RAND_MAX);
		// if (active(v,x) = false) and random <= Pact(u,v)
		if((!v->active[x]) && (rd <= it->second))
		{
			rel(v,R,allNode,Sim,nit,mc,ms); // R <- Rel(v,Sim)

			if(mode == 0)
				UpdatePpref(v,x,R,nit);      // UpdatePpref(v,x,R)

			rd = double(rand())/double(RAND_MAX);
			double Ppref_value = -1.0;
			if(mode == 0)
				Ppref_value = v->Ppref[x];
			else
				Ppref_value = Pminpref;
			if(rd <= Ppref_value) // if random <= Ppref(v,x) then
			{
				v->active[x] = true;   // active(v,x) = true
				union_Inf(newInf,v,x); // newInf <- newInf U {(v,x)}
			}

			EA_pos.clear(); // EA+ <- {}
			EA_neg.clear(); // EA- <- {}
			//for(it2 = allItem.begin();it2 != allItem.end();it2++) // for y in I do
			for(int y=0;y<nit;y++)
			{
				//int y = it2->first;
				double r = R[x][y][0]/(R[x][y][0]+R[x][y][1]);
				rd = double(rand())/double(RAND_MAX);
				if(rd <= r) // if random <= (R ratio) then
				{
					rd = double(rand())/double(RAND_MAX);
					double compare_value = -1.0;
					if(mode == 0)
						compare_value = R[x][y][0];
					else
						compare_value = Pminext * R[x][y][0];

					if(rd <= compare_value) // if random <= R[x][y]['c'] then
					{
						if(mode == 0)
							union_EA(EA_pos,y,'c'); // EA+ <- EA+ U {(y,'c')}
						v->active[y] = true;    // active(v,y) = true
						union_Inf(newInf,v,y);  // newInf <- newInf U {(v,y)}
					}
					else
					{
						if(mode == 0)
							union_EA(EA_neg,y,'c'); // EA- <- EA- U {(y,'c')}
					}
				}
				else
				{
					rd = double(rand())/double(RAND_MAX);
					double compare_value = -1.0;
					if(mode == 0)
						compare_value = R[x][y][1];
					else
						compare_value = Pminext * R[x][y][1];

					if(rd <= compare_value) // if random <= R[x][y]['s'] then
					{
						if(mode == 0)
							union_EA(EA_pos,y,'s'); // EA+ <- EA+ U {(y,'s')}
						v->active[y] = true;    // active(v,y) = true
						union_Inf(newInf,v,y);  // newInf <- newInf U {(v,y)}
					}
					else
					{
						if(mode == 0)
							union_EA(EA_neg,y,'s'); // EA- <- EA- U {(y,'s')}
					}
				}
			}
			// UpdateWeight(v,x,EA+,EA-,Sim)
			if(mode == 0)
				UpdateWeight(v,x,EA_pos,EA_neg,Sim,nit,mc,ms,eta);
		}
	}

	for(int i=0;i<nit;i++)
	{
		for(int j=0;j<nit;j++)
		{
			delete[] R[i][j];
		}
		delete[] R[i];
	}
	delete[] R;
	//return newInf;
}

