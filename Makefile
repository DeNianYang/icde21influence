CFLAGS = -g -Wall
CFLAGS += -O3

all: alg alg_opt

alg: alg.o function.o readfile.o diffusion.o exertInf.o updatePref.o updateWeight.o rel.o greedy.o assignment.o
	g++ -o $@ alg.o function.o readfile.o diffusion.o exertInf.o updatePref.o updateWeight.o rel.o greedy.o assignment.o $(CFLAGS)

alg_opt: alg_opt.o function.o readfile.o diffusion.o exertInf.o updatePref.o updateWeight.o rel.o opt.o
	g++ -pthread -o $@ alg_opt.o function.o readfile.o diffusion.o exertInf.o updatePref.o updateWeight.o rel.o opt.o $(CFLAGS)

alg_opt.o: alg_opt.cpp function.h
	g++ -c alg_opt.cpp $(CFLAGS)

alg.o: alg.cpp function.h
	g++ -c alg.cpp $(CFLAGS)

opt.o: opt.cpp function.h
	g++ -c opt.cpp $(CFLAGS)

assignment.o: assignment.cpp function.h
	g++ -c assignment.cpp $(CFLAGS)

greedy.o: greedy.cpp function.h
	g++ -c greedy.cpp $(CFLAGS)

diffusion.o: diffusion.cpp function.h
	g++ -c diffusion.cpp $(CFLAGS)

exertInf.o: exertInf.cpp function.h
	g++ -c exertInf.cpp $(CFLAGS)

updatePref.o: updatePref.cpp function.h
	g++ -c updatePref.cpp $(CFLAGS)

updateWeight.o: updateWeight.cpp function.h
	g++ -c updateWeight.cpp $(CFLAGS)

rel.o: rel.cpp function.h
	g++ -c rel.cpp $(CFLAGS)

readfile.o: readfile.cpp function.h
	g++ -c readfile.cpp $(CFLAGS)

function.o: function.cpp function.h Data.h
	g++ -c function.cpp $(CFLAGS)

clean:
	rm -rf alg alg_opt *.o
