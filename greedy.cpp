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

void Greedy(vector<Seed*>& GreedyS,map<int,Node*>& allNode,map<int,double>& allItem,double*** Sim,
			int N,int mc,int ms,double eta,double Pminpref,double Pminext,double b,int U_mode,fstream& f)
{
	clock_t time_start = clock();
	// Greedy
	vector<Seed*> U;
	vector<Seed*> S1;
	vector<Seed*> S1_usm;
	vector<Seed*> U_S1;
	vector<Seed*> S2;
	vector<Seed*> S3;
	map<int,Node*>::iterator it;
	map<int,double>::iterator it2;

	int triplet_round = 0;
	for(it = allNode.begin();it != allNode.end();it++)
	{
		for(it2 = allItem.begin();it2 != allItem.end();it2++)
		{
			if(it->second->cost[it2->first] > b) // if cost > b, continue
				continue;
			Seed* newseed;
			if(U_mode == 0) // general case
				newseed = new Seed(it->second,it2->first,1);
			else
			{
				triplet_round++;
				newseed = new Seed(it->second,it2->first,triplet_round);
			}
			U.push_back(newseed);
		}
	}

	//double a = diffusion(allNode,allItem,S,Sim,N,eta,mc,ms,Pminpref,Pminext,mode);
	cout << "   Computing S1 ... First CostGreedy()" << endl;
	CostGreedy(S1,U,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,b);

	if(U_mode == 1)
	{
		cout << "   Computing S1' ... USM();" << endl;
		USM(S1_usm,S1,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext);
	}

	for(size_t i=0;i<U.size();i++)
	{
		bool exist = false;
		for(size_t j=0;j<S1.size();j++)
		{
			if(U[i] == S1[j])
			{
				exist = true;
				break;
			}
		}
		if(!exist)
			U_S1.push_back(U[i]);
	}
	if(U_mode == 1)
	{
		cout << "   Computing S2 ... Second CostGreedy()" << endl;
		CostGreedy(S2,U_S1,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,b);
	}

	cout << "   Computing S3 ... Single Greedy()" << endl;
	SingleGreedy(S3,U,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,b);
	
	/*cout << "S1:" << endl;
	for(size_t i=0;i<S1.size();i++)
		cout << S1[i]->u->num << " " << S1[i]->x << " " << S1[i]->t << endl;
	cout << "S1_usm:" << endl;
	for(size_t i=0;i<S1_usm.size();i++)
		cout << S1_usm[i]->u->num << " " << S1_usm[i]->x << " " << S1_usm[i]->t << endl;
	cout << "S2:" << endl;
	for(size_t i=0;i<S2.size();i++)
		cout << S2[i]->u->num << " " << S2[i]->x << " " << S2[i]->t << endl;
	cout << "S3:" << endl;
	for(size_t i=0;i<S3.size();i++)
		cout << S3[i]->u->num << " " << S3[i]->x << " " << S3[i]->t << endl;*/

	//S <- argmax(sigma(S1),sigma(S1_usm),sigma(S2),sigma(S3));
	max_greedy(GreedyS,S1,S1_usm,S2,S3,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext);
	double sigma_greedy = diffusion(GreedyS,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,0);
	f << "Greedy Seeds" << endl;
	f << sigma_greedy << endl;
	for(size_t i=0;i<GreedyS.size();i++)
		f << GreedyS[i]->u->num << " " << GreedyS[i]->x << " " << GreedyS[i]->t << endl;
	clock_t time_end = clock();
	f << "Time Usage of Greedy: " << double(time_end - time_start)/double(CLOCKS_PER_SEC) << " seconds." << endl;

}

void CostGreedy(vector<Seed*>& S1,vector<Seed*>& U,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext,double b)
{
	S1.clear();
	//vector<Seed*> Su;
	vector<Seed*> S;

	double c = 0.0;

	// while true
	while(true)
	{
		//cout << "nowcost: " << c << " ,budget: " << b << endl;
		double max = 0.0;
		Seed* s = NULL;

		// for (u,x,t) in U do
		for(size_t i=0;i<U.size();i++)
		{
			if(U[i]->u->cost[U[i]->x] > (b-c))
				continue;
			int x = U[i]->x;
			Node* u = U[i]->u;

			bool exist = false;
			for(size_t j=0;j<S.size();j++)
			{
				if(S[j] == U[i])
				{
					exist = true;
					break;
				}
			}
			if(exist)
				continue;

			// val = ( f(S U {(u,x,t)}) - f(S) ) / c_u,x
			double value1 = diffusion(S,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,1);
			S.push_back(U[i]);
			double value2 = diffusion(S,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,1);
			double val = (value2 - value1) / u->cost[x];

			// remove U[i] from S
			diff_Seed(S,U[i]);

			if(val > max)
			{
				max = val;
				s = U[i];
			}
		}
		// if s = {} or c + c_s > b then
		if(s == NULL || (c + (s->u)->cost[s->x] > b))
			break;
		else
		{
			union_Seed(S,s);
			diff_Seed(U,s);
			c = c + (s->u)->cost[s->x];
		}
	}
	
	S1 = S;
}

void USM(vector<Seed*>& S1_usm,vector<Seed*>& S1, map<int,Node*>& allNode,map<int,double>& allItem,
			double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext)
{
	// X <- {}, Y <- S
	vector<Seed*> X;
	vector<Seed*> Y;
	Y = S1;

	// for (u,x,t) in S do
	for(size_t i=0;i<S1.size();i++)
	{
		// x <- max(0,f(X U {(u,x)}) - f(X))
		// y <- max(0,f(Y \ {(u,x)}) - f(Y))
		// p <- 1
		double value1 = diffusion(X,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,1);
		double value2 = diffusion(Y,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,1);
		union_Seed(X,S1[i]);
		diff_Seed(Y,S1[i]);

		double value3 = diffusion(X,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,1);
		double value4 = diffusion(Y,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,1);
		double x = max(0.0,value3 - value1);
		double y = max(0.0,value4 - value2);

		double p = 1.0;
		if(x+y != 0.0) // if x + y != 0 then
			p = x / (x+y);
		double rd = double(rand())/double(RAND_MAX);
		if(rd <= p)
			union_Seed(Y,S1[i]); // X <- XU {(u,x)}
		else
			diff_Seed(X,S1[i]); // Y <- Y - {(u,x)}
	}
	S1_usm = X;
}

void SingleGreedy(vector<Seed*>& S3,vector<Seed*>& U,map<int,Node*>& allNode,map<int,double>& allItem,
					double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext,double b)
{
	double max = 0.0;
	Seed* s = NULL;

	// for (u,x,t) in U and c_u,x <= b
	for(size_t i=0;i<U.size();i++)
	{
		if(U[i]->u->cost[U[i]->x] > b)
			continue;
		vector<Seed*> vec_s;
		vec_s.push_back(U[i]);
		double value = diffusion(vec_s,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,1);
		if(value > max)
		{
			max = value;
			s = U[i];
		}
	}
	S3.clear();
	if(s != NULL)
		S3.push_back(s);
}

void max_greedy(vector<Seed*>& GreedyS,vector<Seed*>& S1,vector<Seed*>& S1_usm,vector<Seed*>& S2,vector<Seed*>& S3,
				map<int,Node*>& allNode,map<int,double>& allItem,double*** Sim,
				int N,int mc,int ms,double eta,double Pminpref,double Pminext)
{
	double value[4];
	value[0] = diffusion(S1,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,0);
	value[1] = diffusion(S1_usm,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,0);
	value[2] = diffusion(S2,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,0);
	value[3] = diffusion(S3,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,0);

	vector<Seed*> max_S;
	int max = -1;
	int max_value = -2147483647;
	string SS[4] = {"S1","S1_usm","S2","S3"};
	for(int i=0;i<4;i++)
	{
		//cout << SS[i] << " value["<<i<<"] = "<<value[i] << endl;
		if(value[i] > max_value)
		{
			max_value = value[i];
			max = i;
		}
	}
	//cout << "The maximum Seeds is ";
	switch(max)
	{
		case 0: /*cout << "S1" << endl;*/ GreedyS = S1; break;
		case 1: /*cout << "S1_usm" << endl;*/ GreedyS = S1_usm; break;
		case 2: /*cout << "S2" << endl;*/ GreedyS = S2; break;
		case 3: /*cout << "S3" << endl;*/ GreedyS = S3; break;
	}
}
