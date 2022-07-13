# icde21influence
This repository consists of the implementation of the paper published in IEEE ICDE 2021: "Influence Maximization Based on Dynamic Personal Perception in Knowledge Graph"

Main Program

1. make clean

2. make all

3. ./alg param.txt => Greedy and Assignment

4. ./alg_opt param_opt.txt => OPT Algorithm

In param.txt
Line 1: <item file path>
Line 2: <graph file path>
Line 3: <meta_c file path>
Line 4: <meta_s file path>
Line 5: <Ppref file path>
Line 6: <Pext file path>
Line 7: <Sim file path>
Line 8: <cost file path>
Line 9: <output seed file path>
Line 10: <N>: an integer
Line 11: <eta>: an double in [0,1]
Line 12: <mc>: meta_c number
Line 13: <ms>: meta_s number
Line 14: <Pminpref>: an double in [0,1]
Line 15: <Pminext>: an double in [0,1]
Line 16: <U_mode>: 1 for counting S1' and S2, 0 otherwise.
Line 17: <b>: budget, an double number
Line 18: <d>: social distance
Line 19: <theta>: an double in [0,1]
Line 20: <T>: an integer


In param_opt.txt
Line 1: <item file path>
Line 2: <graph file path>
Line 3: <meta_c file path>
Line 4: <meta_s file path>
Line 5: <Ppref file path>
Line 6: <Pext file path>
Line 7: <Sim file path>
Line 8: <cost file path>
Line 9: <output seed file path>
Line 10: <N>: an integer
Line 11: <eta>: an double in [0,1]
Line 12: <mc>: meta_c number
Line 13: <ms>: meta_s number
Line 14: <Pminpref>: an double in [0,1]
Line 15: <Pminext>: an double in [0,1]
Line 16: <b>: budget, an double number
Line 17: <NUM_THREAD>: number of thread, an integer

Input Files Format: (all user and item id start from zero.)
1. item file
<item id x>,<weight W_x>

2. graph file (each line means edge u -> v)
<user id u>,<user id v>,<Pact(u,v)>

3.meta_c file (i = 0~(meta_c number)-1)
<user id u>,<item id x>,<i>,W_meta_c(u,x,i)

4.meta_c file (j = 0~(meta_s number)-1)
<user id u>,<item id x>,<j>,W_meta_s(u,x,j)

5.Ppref file
<user id u>,<item id x>,<Ppref(u,x)>

6.Pext file
<user id u>,<item id x>,<item id y>,<Pext(u,x,y)>

7.Sim file (i = 0~meta_c-1 for mc_i, j = meta_c ~ meta_c + meta_s -1 for ms_j)
<item id x>,<item id y>,<i>,<Sim(x,y,mc_i or ms_j)>

8.cost file
<user id u>,<item id x>,<cost(u,x)>
