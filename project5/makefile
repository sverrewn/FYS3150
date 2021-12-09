CC = g++ -std=c++11
FLAGS = -Wextra -Wpedantic -Wall -larmadillo -I include/

main: simulation.out

test: verify_matrix_lib.out

simulation.out: simulation.o matrix_lib.o
	$(CC) $(FLAGS) -o $@ $^

simulation.o: src/simulation.cpp
	$(CC) $(FLAGS) -c $^

verify_matrix_lib.out: verify_matrix_lib.o matrix_lib.o
	$(CC) $(FLAGS) -o $@ $^

verify_matrix_lib.o: src/verify_matrix_lib.cpp
	$(CC) $(FLAGS) -c $^

matrix_lib.o: src/matrix_lib.cpp
	$(CC) $(FLAGS) -c $^