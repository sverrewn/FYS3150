CC = g++ -std=c++11

FLAGS = -Wall -Wextra -Wpedantic -O3

all: poisson_exact.out general_tridiag.out special_tridiag.out

poisson_exact.out: poisson_exact.o algorithms.o
	$(CC) $(FLAGS) -o $@ $^

general_tridiag.out: general_tridiag.o algorithms.o
	$(CC) $(FLAGS) -o $@ $^

special_tridiag.out: special_tridiag.o algorithms.o
	$(CC) $(FLAGS) -o $@ $^

poisson_exact.o: poisson_exact.cpp
	$(CC) $(FLAGS) -c $^

general_tridiag.o: general_tridiag.cpp
	$(CC) $(FLAGS) -c $^

special_tridiag.o: special_tridiag.cpp
	$(CC) $(FLAGS) -c $^

algorithms.o: algorithms.cpp
	$(CC) $(FLAGS) -c $^

clean:
	rm *.out *.o *.bib *.log *.gz *.aux
.PHONY: clean
