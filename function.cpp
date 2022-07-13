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

int str2int(string l)
{
	int ans = 0;
	int digit = l.length();
	for(int i=0;i<digit;i++)
	{
		ans *= 10;
		ans += (l[i]-48);
	}
	return ans;
}

double str2double(string l)
{
	double ans = 0;
	bool digitflag = false;
	int digitcount = 0;
	int len = l.length();
	for(int i=0;i<len;i++)
	{
		if(digitflag)
			digitcount++;
		if(l[i] != '.')
		{
			ans *= 10.0;
			ans += double(l[i]-48);
		}
		else
			digitflag = true;
	}
	for(int i=0;i<digitcount;i++)
		ans /= 10;
	return ans;
}

double get_random()
{
	return double(rand())/double(RAND_MAX);
}

void union_Inf(vector<pair<Node*,int> >& A,Node* v,int x)
{
	for(size_t i=0;i<A.size();i++)
	{
		if(A[i].first == v && A[i].second == x)
			return;
	}
	pair<Node*,int> p(v,x);
	A.push_back(p);
}

void union_EA(vector<pair<int,char> >& A,int y,char c)
{
	for(size_t i=0;i<A.size();i++)
	{
		if(A[i].first == y && A[i].second == c)
			return;
	}
	pair<int,char> p(y,c);
	A.push_back(p);
}

void union_Seed(vector<Seed*>& S,vector<Seed*>& C)
{
	for(size_t i=0;i<C.size();i++)
	{
		bool exist = false;
		for(size_t j=0;j<S.size();j++)
		{
			if(C[i]->u == S[j]->u && C[i]->x == S[j]->x && C[i]->t && S[j]->t)
			{
				exist = true;
				break;
			}
		}
		if(!exist)
			S.push_back(C[i]);
	}
}
void union_Seed(vector<Seed*>& S,Seed* c)
{
	for(size_t i=0;i<S.size();i++)
	{
		if(S[i] == c)
			return;
	}
	S.push_back(c);
}
void union_Node(vector<Node*>& S,vector<Node*>& C)
{
	for(size_t i=0;i<C.size();i++)
	{
		bool exist = false;
		for(size_t j=0;j<S.size();j++)
		{
			if(C[i] == S[j])
			{
				exist = true;
				break;
			}
		}
		if(!exist)
			S.push_back(C[i]);
	}
}
void union_Node(vector<Node*>& S,Node* c)
{
	for(size_t i=0;i<S.size();i++)
	{
		if(S[i] == c)
			return;
	}
	S.push_back(c);
}
void union_vec_Seed(vector<vector<Seed*> >& UC,vector<Seed*>& C)
{
	for(size_t i=0;i<UC.size();i++)
	{
		if(UC[i] == C)
			return;
	}
	UC.push_back(C);
}
void diff_vec_Seed(vector<vector<Seed*> >& UC,vector<Seed*>& C)
{
	int nullcount = 0;
	int notnullcount = 0;
	for(size_t i=0;i<UC.size();i++)
	{
		if(UC[i] == C)
			nullcount++;
		else
		{
			UC[notnullcount] = UC[i];
			notnullcount++;
		}
	}
	for(int i=0;i<nullcount;i++)
		UC.pop_back();
}
void diff_Seed(vector<Seed*>& S,vector<Seed*>& C)
{
	int nullcount = 0;
	int notnullcount = 0;
	for(size_t i=0;i<S.size();i++)
	{
		bool exist = false;
		for(size_t j=0;j<C.size();j++)
		{
			if(S[i]->u == C[j]->u && S[i]->x == C[j]->x && S[i]->t == C[j]->t)
			{
				nullcount++;
				exist = true;
				break;
			}
		}
		if(!exist)
		{
			S[notnullcount] = S[i];
			notnullcount++;
		}
	}
	for(int i=0;i<nullcount;i++)
		S.pop_back();
}
void diff_Seed(vector<Seed*>& S,Seed* c)
{
	int nullcount = 0;
	int notnullcount = 0;
	for(size_t i=0;i<S.size();i++)
	{
		if(S[i] == c)
			nullcount++;
		else
		{
			S[notnullcount] = S[i];
			notnullcount++;
		}
	}
	for(int i=0;i<nullcount;i++)
		S.pop_back();
}
void union_Item(vector<int>& next,int z)
{
	for(size_t i=0;i<next.size();i++)
	{
		if(next[i] == z)
			return;
	}
	next.push_back(z);
}
