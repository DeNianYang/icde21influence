#include<cstdlib>
#include<iostream>
#include<string>
#include<fstream>
#include<queue>
#include<ctime>
#include<vector>
#include<map>
#include<cmath>
#include "function.h"
using namespace std;

void assignment(vector<Seed*>& Smax,vector<Seed*>& S,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int nit,int mc,int ms,double eta,double Pminpref,double Pminext,
				int d,double theta,int T,fstream& f)
{
	clock_t time_start = clock();
	vector<Seed*> S1;
	vector<vector<TM*> > PTMs;
	TMI(S,PTMs,allNode,Sim,nit,mc,ms,d,theta,T); // PTMs <- TMI(G,S,T,Sim)
	vector<Seed*> current_S;
	vector<Seed*> last_S;
	for(size_t i=0;i<PTMs.size();i++) // for TMs in PTMs do
	{
		for(size_t j=0;j<PTMs[i].size();j++) // for TM in TMs
		{
			last_S = S1;
			diff_Seed(last_S,current_S);
			TM* tm = PTMs[i][j];
			int min_T = 1;		// minT <- 1
			if(S1.size() != 0) // if S' != {}
			{
				int temp = 0;
				for(size_t k=0;k<last_S.size();k++) // take the seeds in same TM
				{
					if(last_S[k]->t > temp)
						temp = last_S[k]->t;
				}
				min_T = min(T,temp + 1); // minT <- min{T,max_t{(u,x,t) in S'} + 1}
			}
			int max_T = min(T,min_T + tm->T_TM); // maxT <- min{T,minT + T_TM}

			while(tm->N_TM.size() != 0) // while N_TM != {} do
			{
				int item = DRE(S1,tm,allNode,allItem,Sim,N,nit,mc,ms,eta,Pminpref,Pminext,d); // i <- DRE(TM,Sim,S')
				vector<Seed*> C;
				for(size_t k=0;k<tm->N_TM.size();k++) // C <- {(u,x,t)|x = i and (u,x,t) in N_TM}
				{
					if(tm->N_TM[k]->x == item)
						C.push_back(tm->N_TM[k]);
				}
				diff_Seed(tm->N_TM,C); // N_TM <- N_TM - C
				// TDSI(C,TM,S',T,Sim) modify S'
				last_S = S1;
				diff_Seed(last_S,current_S);
				TDSI(C,S1,last_S,tm,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,T,min_T,max_T);
			}
		}
		current_S = S1;
	}
	clock_t time_end = clock();
	double value1 = diffusion(S,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,0);
	double value2 = diffusion(S1,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,0);
	if(value1 > value2)
		Smax = S;
	else
		Smax = S1;
	//cout << "S' is " << endl;
	//for(size_t i=0;i<S1.size();i++)
	//	cout << S1[i]->u->num << " " << S1[i]->x << " " << S1[i]->t << endl;
	//cout << "S_greedy value = " << value1 << endl;
	//cout << "S' value = " << value2 << endl;
	double sigma_assignment = diffusion(Smax,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext,0);
	f << "Assignment" << endl;
	f << sigma_assignment << endl;
	for(size_t i=0;i<Smax.size();i++)
		f << Smax[i]->u->num << " " << Smax[i]->x << " " << Smax[i]->t << endl;
	
	f << "Time Usage of Assignment: " << double(time_end-time_start)/double(CLOCKS_PER_SEC) << " seconds." << endl;
}

void TMI(vector<Seed*>& S,vector<vector<TM*> >& PTMs,map<int,Node*>& allNode,
			double*** Sim,int nit,int mc,int ms,int d,double theta,int T)
{
	double*** R = new double**[nit];
	for(int i=0;i<nit;i++)
	{
		R[i] = new double*[nit];
		for(int j=0;j<nit;j++)
			R[i][j] = new double[2]();
	}
	
	// R <- RelAvg(G,Sim)
	RelAvg(R,allNode,Sim,nit,mc,ms);
	PTMs.clear(); // PTMs <- {}
	vector<TM*> allTMs; // allTMs <- {}
	vector<vector<Seed*> > UC;
	vector<vector<Seed*> > IC;
	UserClustering(UC,S,allNode,d); // UC <- UserClustering(S)
	for(size_t i=0;i<UC.size();i++) // for uc in UC
	{
		IC.clear();
		ItemClustering(IC,UC[i],R); // IC <- ItemClustering(uc,R)
		for(size_t j=0;j<IC.size();j++) // for ic in IC
		{
			TM* newTM = new TM(IC[j],d); // TM <- InitTM(ic)
			allTMs.push_back(newTM); // allTMs <- allTMs U {TM}
		}
	}

	for(size_t i=0;i<allTMs.size();i++) // for TM1 in allTMs
	{
		for(size_t j=0;j<allTMs.size();j++) // for TM2 in allTMs
		{
			if(i != j) // TM1 != TM2
			{
				// if |VTM1 intersect VTM2| / |VTM1 union VTM2| > theta
				if(get_ratio(allTMs[i],allTMs[j]) > theta)
				{
					// Add an overlap link between TM1 and TM2
					allTMs[i]->link.push_back(allTMs[j]);
					allTMs[j]->link.push_back(allTMs[i]);
				}
			}
		}
	}

	for(size_t i=0;i<allTMs.size();i++)
	{
		if(allTMs[i]->visit)
			continue;
		vector<TM*> CC;
		DFS(CC,allTMs[i]); // for each connected component CC of all TM
		for(size_t j=0;j<CC.size();j++) // for TM1 in CC
		{
			for(size_t k=0;k<CC.size();k++) // for TM2 in CC
			{
				if(j == k) // TM1 != TM2
					continue;
				TM* TM1 = CC[j];
				TM* TM2 = CC[k];
				for(size_t ii=0;ii<TM1->N_TM.size();ii++) // for x in (u,x,t) in N_TM1
				{
					for(size_t jj=0;jj<TM2->N_TM.size();jj++) // for y in (v,y,s) in N_TM2
					{
						int x = TM1->N_TM[ii]->x;
						int y = TM2->N_TM[jj]->x;
						// rTM1 <- rTM1 + R[x][y]['s'] - R[x][y]['c']
						TM1->r_TM = TM1->r_TM + R[x][y][1] - R[x][y][0];
					}
				}
			}
		}
		// Add all TM into TMs in ascending order of rTM
		SortTM(CC);
		int n = 0;
		for(size_t j=0;j<CC.size();j++) // for TM in TMs
			n += int(CC[j]->N_TM.size()); // n <- n + |NTM|
		for(size_t j=0;j<CC.size();j++) // for TM in TMs
			CC[j]->T_TM = floor(double(CC[j]->N_TM.size())*double(T)/double(n)); // TM <- floor(|N_TM|/n * T)
		PTMs.push_back(CC); // PTMs <- PTMs U {TMs}
	}

	// release R
	for(int i=0;i<nit;i++)
	{
		for(int j=0;j<nit;j++)
			delete[] R[i][j];
		delete[] R[i];
	}
	delete[] R;
}

void DFS(vector<TM*>& G,TM* n) // find connected component CC
{
	n->visit = true;
	for(size_t i=0;i<n->link.size();i++)
	{
		if(n->link[i]->visit)
			continue;
		DFS(G,n->link[i]);
	}
	G.push_back(n);
}

void SortTM(vector<TM*>& CC) // sort CC by r_TM
{
	int len = int(CC.size());
	for(int i=0;i<len;i++)
	{
		TM* key = CC[i];
		int j = i - 1;
		while((j > 0) && (key->r_TM < CC[j]->r_TM))
		{
			CC[j+1] = CC[j];
			j--;
		}
		CC[j+1] = key;
	}
}

void UserClustering(vector<vector<Seed*> >& UC,vector<Seed*>& S,map<int,Node*>& allNode,int d)
{
	UC.clear(); // UC <- {}
	for(size_t i=0;i<S.size();i++) // for (u,x,t) in S
	{
		vector<Seed*> C;
		C.push_back(S[i]);     // C <- {(u,x,t)}
		union_vec_Seed(UC,C);  // UC <- UC U {C}
	}
	vector<vector<Seed*> > MC;
	while(UC.size() > 1) // while |UC| > 1
	{
		int min = 2147483647; // min <- inf
		MC.clear(); // MC <- {}
		for(size_t i=0;i<UC.size();i++) // for C1 in UC
		{
			for(size_t j=0;j<UC.size();j++) // for C2 in UC
			{
				if(i != j) // C1 != C2
				{
					int dist = 0; // dist <- 0
					for(size_t k=0;k<UC[i].size();k++) // for u in (u,x,t) in C1
					{
						Node* u = UC[i][k]->u;
						for(size_t m=0;m<UC[j].size();m++) // for v in (v,y,s) in C2
						{
							Node* v = UC[j][m]->u;
							int h = distance(allNode,u,v); // h <- shortest distance between u and v
							if(h > dist)
								dist = h;
						}
					}
					if(dist < min)
					{
						min = dist;
						MC.clear(); // MC <- {C1,C2}
						MC.push_back(UC[i]);
						MC.push_back(UC[j]);
					}
				}
			}
		}
		if(min > d) // if min > d
			break;
		else
		{
			vector<Seed*> Cx; // C' <- {}
			for(size_t i=0;i<MC.size();i++) // for C in MC
			{
				union_Seed(Cx,MC[i]); // C' <- C' U C
				diff_vec_Seed(UC,MC[i]);  // UC <- UC - {C}
			}
			union_vec_Seed(UC,Cx); // UC <- UC U {C'}
		}
	}
}

void ItemClustering(vector<vector<Seed*> >& IC,vector<Seed*>& uc,double*** R)
{
	IC.clear(); // IC <- {}
	for(size_t i=0;i<uc.size();i++) // for (u,x,t) in uc
	{
		vector<Seed*> C;
		C.push_back(uc[i]); // C <- {(u,x,t)}
		union_vec_Seed(IC,C);    // IC <- IC U {C}
	}
	vector<vector<Seed*> > MC;
	while(IC.size() > 1) // while |IC| > 1
	{
		double min = 2147483647.0; // min = inf
		MC.clear(); // MC <- {}
		for(size_t i=0;i<IC.size();i++) // for C1 in IC
		{
			for(size_t j=0;j<IC.size();j++) // for C2 in IC
			{
				if(i != j) // C1 != C2
				{
					double dist = 0.0; // dist <- 0
					for(size_t k=0;k<IC[i].size();k++) // for x in (u,x,t) in C1
					{
						int x = IC[i][k]->x;
						for(size_t m=0;m<IC[j].size();m++) // for y in (v,y,s) in C2
						{
							int y = IC[j][m]->x;
							double r = R[x][y][1] - R[x][y][0]; // h <- R[x][y]['s'] - R[x][y]['c']
							if(r > dist)
								dist = r;
						}
					}
					if(dist < min) // if dist < min
					{
						min = dist;
						MC.clear(); // MC <- {C1,C2}
						MC.push_back(IC[i]);
						MC.push_back(IC[j]);
					}
				}
			}
		}
		if(min > 0) // if min > 0
			break;
		else
		{
			vector<Seed*> Cx; // C' <- {}
			for(size_t i=0;i<MC.size();i++) // for C in MC
			{
				union_Seed(Cx,MC[i]); // C' <- C' U C
				diff_vec_Seed(IC,MC[i]); // IC <- IC - {C}
			}
			union_vec_Seed(IC,Cx); // IC <- IC U {Cx}
		}
	}
}

int DRE(vector<Seed*>& S1,TM* tm,map<int,Node*>& allNode,map<int,double>& allItem,
		double*** Sim,int N,int nit,int mc,int ms,double eta,double Pminpref,double Pminext,int d)
{
	double*** R = new double**[nit];
	for(int i=0;i<nit;i++)
	{
		R[i] = new double*[nit];
		for(int j=0;j<nit;j++)
		{
			R[i][j] = new double[2]();
		}
	}
	//TODO: R <- UpdateRelAvg(V_TM,Sim,S')
	//UpdateRelAvg(tm->V_TM,Sim,S);
	UpdateRelAvg(S1,tm->V_TM,R,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext);
	
	double max_value = -2147483647.0; // max <- -inf
	int item = -1; // i <- {}
	for(size_t i=0;i<tm->N_TM.size();i++) // for x in (u,x,t) in N_TM
	{
		int x = tm->N_TM[i]->x;
		// dr <- ImpactSample(x,R,d,'f') + ImpactSample(x,R,d,'b')
		double dr = ImpactSample(x,R,'f',allItem,N,nit,d) + ImpactSample(x,R,'b',allItem,N,nit,d);
		if(dr >= max_value) // if dr >= max
		{
			max_value = dr;
			item = x;
		}
	}
	// release R
	for(int i=0;i<nit;i++)
	{
		for(int j=0;j<nit;j++)
			delete[] R[i][j];
		delete[] R[i];
	}
	delete[] R;
	
	return item;
}

double ImpactSample(int x,double*** R,char dir,map<int,double>& allItem,int N,int nit,int d)
{
	double val = 0.0; // val <- 0
	vector<int> now;
	vector<int> next;
	vector<int> N_item;
	for(int count=0;count<N;count++) // for i in N
	{
		now.clear();
		now.push_back(x); // now <- {x}
		int h = 0;  // h <- 0
		while(h < d && now.size() != 0) // while h < d and now != {}
		{
			next.clear(); // next <- {}
			for(size_t i=0;i<now.size();i++) // for y in now
			{
				int y = now[i];
				N_item.clear();
				if(dir == 'f') // N <- y's out neighbor
				{
					for(int j=0;j<nit;j++)
					{
						if(R[y][j][0] != 0 || R[y][j][1] != 0)
							N_item.push_back(j);
					}
				}
				else // N <- y's in neighbor
				{
					for(int j=0;j<nit;j++)
					{
						if(R[j][y][0] != 0 || R[j][y][1] != 0)
							N_item.push_back(j);
					}
				}
				for(size_t j=0;j<N_item.size();j++) // for z in N
				{
					int z = N_item[j];
					// assign gain rc rs value
					double gain = (dir == 'f')?(allItem[z]):(1);
					double rc 	= (dir == 'f')?(R[y][z][0]):(R[z][y][0]);
					double rs 	= (dir == 'f')?(R[y][z][1]):(R[z][y][1]);
					double rd 	= double(rand())/double(RAND_MAX);
					if(rd < (rc/(rc+rs)))
						gain = gain * rc;
					else
						gain = gain * (-rs);
					val = val + gain;
					union_Item(next,z); // next <- next U {z}
				}
			}
			now = next;
			h++;
		}
	}
	return (val/double(N));
}

void TDSI(vector<Seed*>& C,vector<Seed*>& S1,vector<Seed*>& last_S,TM* tm,map<int,Node*>& allNode,map<int,double>& allItem,
			double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext,int T,int min_T,int max_T)
{
	int minT = min_T;
	int maxT = max_T;
	if(S1.size() != 0) // if S' != {}
	{
		int temp = 0; // minT <- min{T,max_t{(u,x,t) in S'} + 1}
		for(size_t i=0;i<last_S.size();i++)
		{
			if(last_S[i]->t > temp)
				temp = last_S[i]->t;
		}
		minT = min(T,temp+1);
	}
	pair<double,double> lastIP;
	// (lastI,lastP) <- InfPrefSample(S',Sim,V_TM)
	lastIP = InfPrefSample(S1,tm->V_TM,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext);
	while(C.size() != 0) // while C != {}
	{
		double max = -2147483647; // max <- -Inf
		pair<double,double> maxIP(0.0,0.0); // (maxI,maxP) <- (0,0)
		Seed* c = NULL; // c <- {}
		Seed* s = NULL; // s <- {}
		for(size_t i=0;i<C.size();i++) // for (u,x,t) in C
		{
			Node* u = C[i]->u;
			int x = C[i]->x;
			if(minT >=2) // if minT >= 2
			{
				// (i,p) <- InfPrefSample(S' U {(u,x,minT-1)}, Sim)
				pair<double,double> IP;
				Seed* newseed = new Seed(u,x,minT-1);
				S1.push_back(newseed);
				IP = InfPrefSample(S1,tm->V_TM,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext);
				double SI1 = (IP.first - lastIP.first) + (double(T-minT)/double(T)) * (IP.second - lastIP.second);
				diff_Seed(S1,newseed); // S' <- S' - {(u,x,t)}
				if(SI1 > max) // if SI1 > max
				{
					max = SI1;
					maxIP = IP;
					c = C[i];
					s = newseed;
				}
			}
			if(minT <= maxT) // if minT <= maxT
			{
				// (i,p) <- InfPrefSample(S' U {(u,x,minT)}, Sim)
				pair<double,double> IP;
				Seed* newseed = new Seed(u,x,minT);
				S1.push_back(newseed);
				IP = InfPrefSample(S1,tm->V_TM,allNode,allItem,Sim,N,mc,ms,eta,Pminpref,Pminext);
				double SI2 = (IP.first - lastIP.first) + (double(T-minT+1)/double(T)) * (IP.second - lastIP.second);
				diff_Seed(S1,newseed); // S' <- S' - newseed
				if(SI2 > max) // if SI2 > max
				{
					max = SI2;
					maxIP = IP;
					c = C[i];
					s = newseed;
				}
			}
		}
		union_Seed(S1,s); // S' <- S' U {s}
		diff_Seed(C,c); // C <- C - {c}
		lastIP = maxIP; // (lastI,lastP) <- (maxI,maxP)
		int temp = 0;
		for(size_t i=0;i<S1.size();i++)
		{
			if(S1[i]->t > temp)
				temp = S1[i]->t;
		}
		minT = min(T,temp + 1); // minT <- min{max_t{(u,x,t) in S'} + 1}
	}
}

pair<double,double> InfPrefSample(vector<Seed*>& Seeds,vector<Node*>& vtm,map<int,Node*>& allNode,map<int,double>& allItem,
				double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext)
{
	vector<vector<pair<Node*,int> > > S;
	convert_seed(Seeds,S);

	double sigma_S = 0.0; // sigma(S) <- 0
	double pref = 0.0;    // pref <- 0
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
					exertInf(u,x,newInf,allNode,allItem,Sim,mc,ms,eta,Pminpref,Pminext,0);
					for(size_t j=0;j<newInf.size();j++) // for (v,y) in newInf do
					{
						Node* v = newInf[j].first;
						int y = newInf[j].second;
						bool exist = false; // if v in V_TM do
						for(size_t k=0;k<vtm.size();k++)
						{
							if(v == vtm[k])
							{
								exist = true;
								break;
							}
						}
						if(exist)
							sigma_S += allItem[y]; // sigma_S <- sigma_S + W_y
						else
						{
							// Remove (v,y) from Inf
							newInf[j] = newInf[newInf.size()-1];
							newInf.pop_back();
							j--;
						}
					}
					
					// nextInf <- nextInf U newInf
					for(size_t j=0;j<newInf.size();j++)
						union_Inf(nextInf,newInf[j].first,newInf[j].second);
				}
				Inf = nextInf; // Inf <- nextInf
			}
			//cout << "sigma_S = " << sigma_S << endl;
		}
		for(size_t i=0;i<vtm.size();i++) // for u in V_TM 
		{
			// R <- Rel(u,Sim)
			double*** R = new double**[nit];
			for(int j=0;j<nit;j++)
			{
				R[j] = new double*[nit];
				for(int k=0;k<nit;k++)
					R[j][k] = new double[2];
			}

			Node* u = vtm[i];
			rel(u,R,allNode,Sim,nit,mc,ms);
			map<int,double>::iterator it2;
			for(it2 = allItem.begin();it2 != allItem.end();it2++) // for x in I do
			{
				int x = it2->first;
				if(!u->active[x]) //  if active(u,x) != true do
				{
					UpdatePpref(u,x,R,nit);      // UpdatePpref(u,x,R)
					pref = pref + u->Ppref[x] * allItem[x]; // pref <- pref + Ppref(u,x) * W_x
				}
			}

			// release R
			for(int j=0;j<nit;j++)
			{
				for(int k=0;k<nit;k++)
					delete[] R[j][k];
				delete[] R[j];
			}
			delete[] R;
		}

		//cout << "new_sigma_S = " << sigma_S - now_sigma_S << endl;
		// Reset all
		for(it = allNode.begin();it != allNode.end();it++)
			it->second->reset(nit,mc,ms,0);
	}
	pair<double,double> pa(sigma_S/double(N),pref/double(N));
	return pa;
}

int distance(map<int,Node*>& allNode,Node* u,Node* v)
{
	distance(u);
	int d = v->dist;
	map<int,Node*>::iterator it;
	for(it = allNode.begin();it != allNode.end();it++)
		it->second->dist = -1;
	return d;
}
void distance(Node* u)
{
	u->dist = 0;
	queue<Node*> q;
	q.push(u);
	map<Node*,double>::iterator it;
	while(q.size() != 0)
	{
		Node* a = q.front();
		q.pop();
		for(it = a->nout.begin();it != a->nout.end();it++)
		{
			Node* b = it->first;
			if(b->dist < 0)
			{
				b->dist = a->dist + 1;
				q.push(b);
			}
		}
	}
}

void RelAvg(double*** R,map<int,Node*>& allNode,double*** Sim,int nit,int mc,int ms)
{
	for(int i=0;i<nit;i++)
	{
		for(int j=0;j<nit;j++) // for x,y in I do
		{
			R[i][j][0] = 0;
			R[i][j][1] = 0;
			map<int,Node*>::iterator it;
			for(it = allNode.begin();it != allNode.end();it++)
			{
				Node* v = it->second;
				for(int k=0;k<mc;k++) // for meta_graph complementary rel. do
					R[i][j][0] += (v->W_meta_c[k] * Sim[i][j][k]);
				for(int k=0;k<ms;k++) // for meta_graph substitute rel. do
					R[i][j][1] += (v->W_meta_s[k] * Sim[i][j][mc+k]);
			}
			R[i][j][0] = R[i][j][0] / double(allNode.size());
			R[i][j][1] = R[i][j][1] / double(allNode.size());
		}
	}
}

double get_ratio(TM* a,TM* b)
{
	vector<Node*> union_set;
	vector<Node*> intersect_set;
	for(size_t i=0;i<a->V_TM.size();i++)
		union_Node(union_set,a->V_TM[i]);
	for(size_t i=0;i<b->V_TM.size();i++)
		union_Node(union_set,b->V_TM[i]);
	for(size_t i=0;i<a->V_TM.size();i++)
	{
		for(size_t j=0;j<b->V_TM.size();j++)
		{
			if(a->V_TM[i] == b->V_TM[j])
				union_Node(intersect_set,a->V_TM[i]);
		}
	}
	return double(intersect_set.size())/double(union_set.size());
}
void UpdateRelAvg(vector<Seed*>& Seeds,vector<Node*>& VTM,double*** R,map<int,Node*>& allNode,map<int,double>& allItem,
					double*** Sim,int N,int mc,int ms,double eta,double Pminpref,double Pminext)
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

	for(int i=0;i<nit;i++)
	{
		for(int j=0;j<nit;j++)
		{
			R[i][j][0] = 0.0;
			R[i][j][1] = 0.0;
		}
	}

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
					exertInf(u,x,newInf,allNode,allItem,Sim,mc,ms,eta,Pminpref,Pminext,0);
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
		//compute R

		for(int i=0;i<nit;i++)
		{
			for(int j=0;j<nit;j++) // for x,y in I do
			{
				for(size_t t=0;t<VTM.size();t++)
				{
					Node* v = VTM[t];
					for(int k=0;k<mc;k++) // for meta_graph complementary rel. do
						R[i][j][0] += (v->W_meta_c[k] * Sim[i][j][k]);
					for(int k=0;k<ms;k++) // for meta_graph substitute rel. do
						R[i][j][1] += (v->W_meta_s[k] * Sim[i][j][mc+k]);
				}
			}
		}


		// Reset all
		for(it = allNode.begin();it != allNode.end();it++)
			it->second->reset(nit,mc,ms,0);
	}
	
	for(int i=0;i<nit;i++)
	{
		for(int j=0;j<nit;j++)
		{
			R[i][j][0] = R[i][j][0]/double(VTM.size())/double(N);
			R[i][j][1] = R[i][j][1]/double(VTM.size())/double(N);
		}
	}

}
