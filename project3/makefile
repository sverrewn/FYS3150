CC = g++ -std=c++11
FLAGS = -O3 -larmadillo -I include/


run_simulation.out: run_simulation.o penning_trap.o particle.o
	$(CC) $(FLAGS) -o $@ $^

test_penning.out: test_penning.o penning_trap.o particle.o
	$(CC) $(FLAGS) -o $@ $^

run_simulation.o: src/run_simulation.cpp
	$(CC) $(FLAGS) -c $^

test_penning.o: src/test_penning.cpp
	$(CC) $(FLAGS) -c $^

penning_trap.o: src/penning_trap.cpp
	$(CC) $(FLAGS) -c $^

particle.o: src/particle.cpp
	$(CC) $(FLAGS) -c $^

clean:
	rm *.o
.PHONY: clean
