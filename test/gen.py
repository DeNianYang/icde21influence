from random import randint

I = 10
U = 100
mc = 5
ms = 5
seedN = 15
seedNN = 3
p = 1
q = 150
r = 1000
graph_dense = 0.3

f = open("item.txt","w")
for i in range(I):
	f.write(str(i)+','+str(randint(p,q)/r)+'\n')
f.close()

f = open("graph.txt",'w')
for i in range(U):
	for j in range(U):
		if i != j:
			if randint(1,100)/100 < graph_dense:
				f.write(str(i)+','+str(j)+','+str(randint(p,q)/r)+'\n')
f.close()

f = open("seed.txt",'w')
used = []
for i in range(seedN):
	for j in range(seedNN):
		a = randint(0,U-1)
		b = randint(0,I-1)
		if (a,b) not in used:
			used.append((a,b))
			f.write(str(a)+','+str(b)+','+str(i)+'\n')
f.close()

f = open("mc.txt",'w')
for i in range(U):
	for j in range(mc):
		f.write(str(i)+','+str(j)+','+str(randint(p,q)/r)+'\n')
f.close()

f = open("ms.txt",'w')
for i in range(U):
	for j in range(ms):
		f.write(str(i)+','+str(j)+','+str(randint(p,q)/r)+'\n')
f.close()

f = open("pref.txt",'w')
for i in range(U):
	for j in range(I):
		f.write(str(i)+','+str(j)+','+str(randint(p,q)/r)+'\n')
f.close()

f = open("pext.txt",'w')
for i in range(U):
	for j in range(I):
		for k in range(I):
			f.write(str(i)+','+str(j)+','+str(k)+','+str(randint(p,q)/r)+'\n')
f.close()

f = open("sim.txt","w")
for i in range(I):
	for j in range(I):
		for k in range(mc+ms):
			f.write(str(i)+','+str(j)+','+str(k)+','+str(randint(p,q)/r)+'\n')
f.close()

Cost = [0.2,0.3,0.4,0.5,0.6]

f = open("cost.txt",'w')
for i in range(U):
	for j in range(I):
		f.write(str(i)+','+str(j)+','+str(Cost[randint(1,10000)%5])+'\n')
