CC = g++ -std=c++11
FLAGS = -O3 -larmadillo -I include/


all: test_ising.out

test_ising.out: test_ising.o lattice.o
	$(CC) $(FLAGS) -fopenmp -o $@ $^

test_ising.o: src/test_ising.cpp
	$(CC) $(FLAGS) -fopenmp -c $^

lattice.o: src/lattice.cpp
	$(CC) $(FLAGS) -c $^


clean:
	rm *.o
.PHONY: clean
