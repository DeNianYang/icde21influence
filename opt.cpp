#include<cstdlib>
#include<iostream>
#include<string>
#include<fstream>
#include<queue>
#include<ctime>
#include<vector>
#include<map>
#include<thread>
#include<mutex>
#include "function.h"
using namespace std;

mutex mu;

void test_thread(void* p)
{
	Param* param = (Param*)p;
	vector<Seed*> &U = *(param->U);
	map<int,Node*> &allNode = *(param->allNode);
	map<int,double> &allItem = *(param->allItem);
	vector<Seed*> &S1 = *(param->S1);
	vector<Seed*> &S = *(param->S);
	vector<vector<Seed*> > &now = *(param->now);
	double*** Sim = (param->Sim);
	int N = (param->N);
	int mc = (param->mc);
	int ms = (param->ms);
	double eta = (param->eta);
	int T = (param->T);
	double b = (param->b);
	double &max = *(param->max);
	vector<Seed*> S2;
	vector<vector<Seed*> > next;
	for(size_t j=0;j<U.size();j++) // for e = (u,x) in U do
	{
		//cout << "    " << j+1 << "/" << U.size() << endl;
		Node* u = U[j]->u;
		int x = U[j]->x;
		double cost_value = u->cost[x];
		for(size_t k=0;k<S1.size();k++)
			cost_value += (S1[k]->u)->cost[S1[k]->x];
		if(cost_value <= b) // if c_u,x + sum_(v,y)in_S' {c_v,y} <=b then
		{
			next.clear();
			EnumerateRounds(next,S1,u,x,T); // next <- EnumerateRounds(S',e,T)
			for(size_t l=0;l<next.size();l++) // for S'' in next do
			{
				S2 = next[l];
				double value = diffusion(S2,allNode,allItem,Sim,N,mc,ms,eta,0.0,0.0,0);
				mu.lock();
				if(value > max) // if sigma(S'') > max then
				{
					max = value; // max <- sigma(S'')
					S = S2; // S <- S''
				}
				union_vec_Seed(now,S2); // now <- now U {S''}
				mu.unlock();
			}
		}
	}
	cout << "Process No." << param->proc_num << " finish!" << endl;
}

void opt(vector<Seed*>& S,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int mc,int ms,double eta,double b,int T,fstream& f,int NUM_THREAD)
{
	cout << "OPT ... " << endl;
	clock_t time_start = clock();
	map<int,Node*>::iterator it1;
	map<int,double>::iterator it2;

	// U <- {(u,x)| u in allUser, x in allItem, c_u,x <= b}
	vector<Seed*> U;
	for(it1 = allNode.begin();it1 != allNode.end();it1++)
	{
		for(it2 = allItem.begin();it2 != allItem.end();it2++)
		{
			if(it1->second->cost[it2->first] <= b)
			{
				Seed* newseed = new Seed(it1->second,it2->first,-1);
				U.push_back(newseed);
			}
		}
	}
	cout << "U.size() = " << U.size() << endl;
	double max = 0.0; // max <- 0
	S.clear(); // S <- {}
	vector<vector<Seed*> > now; // now <- {}
	vector<Seed*> S1;
	vector<Seed*> S2;
	for(size_t i=0;i<U.size();i++) // for e = (u,x) in U do
	{
		Seed* newseed = new Seed(U[i]->u,U[i]->x,1);
		S1.clear();
		S1.push_back(newseed); // S' <- {(u,x,1)}
		double value = diffusion(S1,allNode,allItem,Sim,N,mc,ms,eta,0.0,0.0,0);
		if(value > max) // if sigma(S') > max then
		{
			max = value; // max <- sigma(S')
			S = S1; // S <- S'
		}
		union_vec_Seed(now,S1); // now <- now U {S'}
	}

	cout << "Time Usage: " << double(clock() - time_start)/double(CLOCKS_PER_SEC) << " seconds." << endl;

	vector<vector<Seed*> > next;
	while(now.size() != 0) // while now != {} do
	{
		cout << "now.size() = " << now.size() << endl;
		for(size_t i=0;i<now.size();i++) // for S' in now do
		{
			cout << i+1 << "/" << now.size() << endl;
			S1 = now[i];

			// using thread
			vector<vector<Seed*> > partial_U;
			int num_thread = NUM_THREAD;
			for(int j=0;j<num_thread;j++)
			{
				vector<Seed*> temp;
				partial_U.push_back(temp);
			}
			int count = 0;
			int now_count = 0;
			for(size_t j=0;j<U.size();j++)
			{
				if(count > double(U.size())/double(num_thread))
				{
					now_count++;
					count = 0;
				}
				partial_U[now_count].push_back(U[j]);
				count++;
			}

			vector<thread> vec_thread;
			for(int j=0;j<num_thread;j++)
			{
				Param* param = new Param();
				param->U = &partial_U[j];
				param->allNode = &allNode;
				param->allItem = &allItem;
				param->S1 = &S1;
				param->S = &S;
				param->now = &now;
				param->Sim = Sim;
				param->N = N;
				param->mc = mc;
				param->ms = ms;
				param->eta = eta;
				param->T = T;
				param->b = b;
				param->max = &max;
				param->proc_num = j+1;
				vec_thread.push_back(thread(&test_thread,param));
			}
			for(int j=0;j<num_thread;j++)
				vec_thread[j].join();


			/*for(size_t j=0;j<U.size();j++) // for e = (u,x) in U do
			{
				//cout << "    " << j+1 << "/" << U.size() << endl;
				Node* u = U[j]->u;
				int x = U[j]->x;
				double cost_value = u->cost[x];
				for(size_t k=0;k<S1.size();k++)
					cost_value += (S1[k]->u)->cost[S1[k]->x];
				if(cost_value <= b) // if c_u,x + sum_(v,y)in_S' {c_v,y} <=b then
				{
					next.clear();
					EnumerateRounds(next,S1,u,x,T); // next <- EnumerateRounds(S',e,T)
					for(size_t l=0;l<next.size();l++) // for S'' in next do
					{
						S2 = next[l];
						double value = diffusion(S2,allNode,allItem,Sim,N,mc,ms,eta,0.0,0.0,0);
						if(value > max) // if sigma(S'') > max then
						{
							max = value; // max <- sigma(S'')
							S = S2; // S <- S''
						}
						union_vec_Seed(now,S2); // now <- now U {S''}
					}
				}

			}*/
			diff_vec_Seed(now,S1);
			i--;
		}
	}
	double sigma_S = diffusion(S,allNode,allItem,Sim,N,mc,ms,eta,0.0,0.0,0);
	clock_t time_end = clock();
	f << "OPT" << endl;
	f << sigma_S << endl;
	for(size_t i=0;i<S.size();i++)
		f << S[i]->u->num << " " << S[i]->x << " " << S[i]->t << endl;
	f << "Time Usage of OPT: " << double(time_end-time_start)/double(CLOCKS_PER_SEC) << " seconds." << endl;
}

void EnumerateRounds(vector<vector<Seed*> >& next,vector<Seed*>& S1,Node* u,int x,int T)
{
	next.clear(); // next <- {}
	int maxT = -2147483647;
	for(size_t i=0;i<S1.size();i++) // maxT = max_t{t|(u,x,t') in S'}
	{
		if(S1[i]->t > maxT)
			maxT = S1[i]->t;
	}
	bool ins = false; // if maxT <T then ins = true ,else ins = false
	if(maxT < T)
		ins = true;

	vector<Seed*> S2; // S''
	vector<Seed*> S11; // S1
	vector<Seed*> S22; // S2
	for(int t=1;t<=maxT;t++) // for t = 1,2, ... , maxT
	{
		S2.clear();
		S11.clear();
		S22.clear();
		if(!seed_in_vec(u,x,t,S1)) // if (u,x,t) not in S'
		{
			S2 = S1; // S'' <- S' U {(u,x,y)}
			Seed* newseed = new Seed(u,x,t);
			union_Seed(S2,newseed);
			union_vec_Seed(next,S2); // next <- next U {S''}
		}
		if(ins) // if ins then
		{
			for(size_t i=0;i<S1.size();i++) // S1 <- {(u',x',t')|(u',x',t') in S',t' < t}
			{
				if(S1[i]->t < t)
					S11.push_back(S1[i]);
			}
			S2 = S11; // S'' <- S1 U {(u,x,t)}
			Seed* newseed = new Seed(u,x,t);
			union_Seed(S2,newseed);
			S22 = S1; // S2 <- S' - S1
			diff_Seed(S22,S11);
			for(size_t i=0;i<S22.size();i++) // for (u',x',t') in S2 do 
			{
				Seed* newseed2 = new Seed(S22[i]->u,S22[i]->x,S22[i]->t + 1); // S'' <- S'' U {(u',x',t'+1)}
				union_Seed(S2,newseed2);
			}
			union_vec_Seed(next,S2); // next <- next U {S''}
		}
	}
	if(ins) // if ins then
	{
		S2 = S1; // S'' <- S' U (u,x,maxT + 1)
		Seed* newseed = new Seed(u,x,maxT+1);
		union_Seed(S2,newseed);
		union_vec_Seed(next,S2); // next <- next U {S''}
	}
}

bool seed_in_vec(Node* u,int x,int t,vector<Seed*>& S1)
{
	for(size_t i=0;i<S1.size();i++)
	{
		if(S1[i]->u == u && S1[i]->x == x && S1[i]->t == t)
			return true;
	}
	return false;
}
