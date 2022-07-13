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

void UpdateWeight(Node* v,int x,vector<pair<int,char> >& EA_pos,vector<pair<int,char> >& EA_neg,
					double*** Sim,int nit,int mc,int ms,double eta)
{
	double W_c = 0.0; // W_c = 0
	for(int m=0;m<mc;m++) // for meta_graph complemetary rel. do
	{
		// Yic <- {y|Sim(x,y,mci) > 0}
		vector<int> Y_i_c;
		for(int i=0;i<nit;i++)
		{
			if(Sim[x][i][m] > 0)
				Y_i_c.push_back(i);
		}

		// nc+ <- 0, nc- <- 0
		int nc_pos = 0;
		int nc_neg = 0;
		for(size_t i=0;i<EA_pos.size();i++) // for (y,'c') in EA+ do
		{
			if(EA_pos[i].second == 'c')     // only for 'c'
			{
				int y = EA_pos[i].first;
				if(exist_item(Y_i_c,y)) // if y in Yic
					nc_pos++;			// nc+ <- nc+ + 1
			}
		}
		for(size_t i=0;i<EA_neg.size();i++) // for (y,'c') in EA- do
		{
			if(EA_neg[i].second == 'c')     // only for 'c'
			{
				int y = EA_neg[i].first;
				if(exist_item(Y_i_c,y)) // if y in Yic
					nc_neg++;           // nc- <- nc- + 1
			}
		}
		// W_meta(v,mci) <- max{0,W_meta(v,mci) + eta(nc+ - nc-)}
		v->W_meta_c[m] = max(0.0,v->W_meta_c[m] + eta*(nc_pos - nc_neg));
		W_c += v->W_meta_c[m]; // W_c += W_meta(v,mci)
	}
	for(int m=0;m<mc;m++) // for meta_graph for complementray rel. do
		v->W_meta_c[m] = v->W_meta_c[m] / W_c; // W_meta(v,mci) <- W_meta(v,mci) / W_c

	double W_s = 0.0; // W_s = 0
	for(int m=0;m<ms;m++) // for meta_graph for substitue rel. do
	{
		// Yjs <- {y|Sim(x,y,msj) > 0}
		vector<int> Y_j_s;
		for(int i=0;i<nit;i++)
		{
			if(Sim[x][i][mc+m] > 0)
				Y_j_s.push_back(i);
		}

		// ns+ <- 0, ns- <- 0
		int ns_pos = 0;
		int ns_neg = 0;
		for(size_t i=0;i<EA_pos.size();i++) // for (y,'s') in EA+ do
		{
			if(EA_pos[i].second == 's')     // only for 's'
			{
				int y = EA_pos[i].first;
				if(exist_item(Y_j_s,y)) // if y in Yjs
					ns_pos++;           // ns+ <- ns+ + 1
			}
		}
		for(size_t i=0;i<EA_neg.size();i++) // for (y,'s') in EA- do
		{
			if(EA_neg[i].second == 's')     // only for 's'
			{
				int y = EA_neg[i].first;
				if(exist_item(Y_j_s,y)) // if y in Yjs
					ns_neg++;			// ns- <- ns- + 1
			}
		}
		// W_meta(v,msj) <- max{0,W_meta(v,msj) + eta*(ns+ - ns-)}
		v->W_meta_s[m] = max(0.0,v->W_meta_s[m] + eta*(ns_pos - ns_neg));
		W_s += v->W_meta_s[m]; // W_s <- W_s + W_meta(v,msj)
	}
	for(int m=0;m<ms;m++) // for meta_graph for substitute rel. do
		v->W_meta_s[m] = v->W_meta_s[m] / W_s; // W_meta(v,msj) <- W_meta(v,msj) / W_s
}

bool exist_item(vector<int>& A,int a)
{
	for(size_t i=0;i<A.size();i++)
	{
		if(A[i] == a)
			return true;
	}
	return false;
}
