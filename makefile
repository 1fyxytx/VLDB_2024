CCOMPILE=mpicxx
PLATFORM=Linux-amd64-64
CPPFLAGS=-I$(HOME)/blogel

all: run

run: run.cpp
	$(CCOMPILE) -O2 -std=c++11 -fopenmp run.cpp $(CPPFLAGS) -o run 

clean:
	-rm run