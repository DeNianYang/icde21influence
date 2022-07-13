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

double diffusion(vector<Seed*>& Seeds,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext,int mode)
{
	vector<vector<pair<Node*,int> > > S;
	convert_seed(Seeds,S);

	double sigma_S = 0.0; // sigma(S) <- 0
	int nit = int(allItem.size());
	vector<pair<Node*,int> > Inf;
	vector<pair<Node*,int> > nextInf;
	vector<pair<Node*,int> > newInf;
	map<int,Node*>::iterator it;
	Node* u;
	int x;

	for(int p=0;p<N;p++)  // for n in N do
	{
		//double now_sigma_S = sigma_S;
		for(size_t t=0;t<S.size();t++) // for t in T do
		{
			Inf = S[t]; // Inf <- St
			for(size_t i=0;i<Inf.size();i++)
			{
				u = Inf[i].first;
				x = Inf[i].second;
				u->active[x] = true;
			}
			while(Inf.size() != 0)  // while Inf is not {}
			{
				nextInf.clear(); // nextInf <- {}
				for(size_t i=0;i<Inf.size();i++)  // for (u,x) in Inf do
				{
					// newInf <- ExertInf((u,x))
					u = Inf[i].first;
					x = Inf[i].second;
					newInf.clear();
					exertInf(u,x,newInf,allNode,allItem,Sim,mc,ms,eta,Pminpref,Pminext,mode);
					for(size_t j=0;j<newInf.size();j++) // for (v,y) in newInf do
						sigma_S += allItem[newInf[j].second];

					// nextInf <- nextInf U newInf
					for(size_t j=0;j<newInf.size();j++)
						union_Inf(nextInf,newInf[j].first,newInf[j].second);
				}
				Inf = nextInf; // Inf <- nextInf
			}
			//cout << "sigma_S = " << sigma_S << endl;
		}
		//cout << "new_sigma_S = " << sigma_S - now_sigma_S << endl;
		// Reset all
		for(it = allNode.begin();it != allNode.end();it++)
			it->second->reset(nit,mc,ms,mode);
	}

	return sigma_S/double(N);
}

void convert_seed(vector<Seed*>& Seed,vector<vector<pair<Node*,int> > >& S)
{
	//InsertionSort(Seed,int(Seed.size()));
	QuickSort(Seed,0,int(Seed.size()-1));
	int max_t = -1;
	vector<pair<Node*,int> > seed_t;
	for(size_t i=0;i<Seed.size();i++)
	{
		if(Seed[i]->t > max_t)
		{
			if(seed_t.size() > 0)
				S.push_back(seed_t);
			max_t = Seed[i]->t;
			seed_t.clear();
		}
		pair<Node*,int> p(Seed[i]->u,Seed[i]->x);
		seed_t.push_back(p);
	}
	S.push_back(seed_t);
}
void InsertionSort(vector<Seed*>& A,int len)
{
	for(int i=1;i<len;i++)
	{
		Seed* key = A[i];
		int j = i-1;
		while((j >= 0) && (key->t < A[j]->t))
		{
			A[j+1] = A[j];
			j--;
		}
		A[j+1] = key;
	}
}
void QuickSort(vector<Seed*>& A,int front,int end)
{
	if(front < end)
	{
		int pivot = Partition(A,front,end);
		QuickSort(A,front,pivot-1);
		QuickSort(A,pivot+1,end);
	}
}
int Partition(vector<Seed*>& A,int front,int end)
{
	int pivot = A[end]->t;
	int i = front - 1;
	for(int j=front;j<end;j++)
	{
		if(A[j]->t < pivot)
		{
			i++;
			Seed* temp = A[i];
			A[i] = A[j];
			A[j] = temp;
		}
	}
	i++;
	Seed* temp = A[i];
	A[i] = A[end];
	A[end] = temp;
	return i;
}
