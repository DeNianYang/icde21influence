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

void UpdatePpref(Node* v,int x,double*** R,int nit)
{
	// max_c = 0, max_s = 0
	double max_c = 0;
	double max_s = 0;

	//for(size_t i=0;i<IA_v.size();i++) // for y in IA(v) do
	for(int y=0;y<nit;y++)
	{
		if(!v->active[y]) // only for active(v,y) == true
			continue;
		double r = R[x][y][0]/(R[x][y][0]+R[x][y][1]);
		double rd = double(rand())/double(RAND_MAX);
		if(rd <= r) // if random <= (R ratio) then
		{
			if(R[x][y][0] > max_c)  // if R[x][y]['c'] > max_c then
				max_c = R[x][y][0]; // max_c <- R[x][y]['c']
		}
		else
		{
			if(R[x][y][1] > max_s)  // if R[x][y]['s'] > max_s then
				max_s = R[x][y][1]; // max_s <- R[x][y]['s']
		}
	}
	if(max_c > max_s)						  // if max_c > max_s
		v->Ppref[x] = max(v->Ppref[x],max_c);   // Ppref(v,x) <- max{Ppref(v,x),max_c}
	else
		v->Ppref[x] = min(v->Ppref[x],1-max_s); // Ppref(v,x) <- min{Ppref(v,x),1-max_s}
}

