CC = g++ -std=c++11

WARN = -Wall -Wextra -Wpedantic

test: poisson_exact.test general_tridiag.test

release: poisson_exact.out general_tridiag.out

poisson_exact.out: poisson_exact.cpp
	$(CC) -Ofast -o $@ $^

general_tridiag.out: general_tridiag.cpp
	$(CC) -Ofast -o $@ $^

general_tridiag.test: general_tridiag.cpp
	$(CC) $(WARN) -o $@ $^

poisson_exact.test: poisson_exact.cpp
	$(CC) $(WARN) -o $@ $^

clean:
	rm *.out *.test *.dat *.bib *.log *.gz *.aux
.PHONY: clean