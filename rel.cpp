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

void rel(Node* v,double*** R,map<int,Node*>& allNode,double*** Sim,int nit,int mc,int ms)
{
	for(int i=0;i<nit;i++)
	{
		for(int j=0;j<nit;j++) // for x,y in I do
		{
			R[i][j][0] = 0;
			for(int k=0;k<mc;k++) // for meta_graph complementary rel. do
				R[i][j][0] += (v->W_meta_c[k] * Sim[i][j][k]); // R[x][y]['c'] += W_meta(v,mci)*Sim(x,y,mci)
			R[i][j][1] = 0;
			for(int k=0;k<ms;k++) // for meta_graph substitute rel. do
				R[i][j][1] += (v->W_meta_s[k] * Sim[i][j][mc+k]); // R[x][y]['s'] += W_meta(v,msj)*Sim(x,y,msj)
		}
	}
}

