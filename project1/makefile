CC = g++ -std=c++11

WARN = -Wall -Wextra -Wpedantic

all: poisson_exact.out

poisson_exact.out: poisson_exact.cpp
	$(CC) $(WARN) -o $@ $^

clean:
	rm *.out *.dat
.PHONY: clean