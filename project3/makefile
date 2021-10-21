CC = g++ -std=c++11
FLAGS = -Wall -O3 -Wextra -Wpedantic -larmadillo -I include/

test_penning.out: test_penning.o penning_trap.o particle.o
	$(CC) $(FLAGS) -o $@ $^
	rm *.o

test_penning.o: src/test_penning.cpp
	$(CC) $(FLAGS) -c $^

penning_trap.o: src/penning_trap.cpp
	$(CC) $(FLAGS) -c $^

particle.o: src/particle.cpp
	$(CC) $(FLAGS) -c $^

clean:
	rm *.out
.PHONY: clean