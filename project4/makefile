CC = g++ -std=c++11
FLAGS = -O3 -larmadillo -I include/


all: run_ising.out

run_ising.out: run_ising.o lattice.o
	$(CC) $(FLAGS) -fopenmp -o $@ $^

run_ising.o: src/run_ising.cpp
	$(CC) $(FLAGS) -fopenmp -c $^

lattice.o: src/lattice.cpp
	$(CC) $(FLAGS) -c $^


clean:
	rm *.o
.PHONY: clean
