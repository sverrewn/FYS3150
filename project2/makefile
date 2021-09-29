CC = g++ -std=c++11

FLAGS = -Wall -Wextra -Wpedantic -O2 -larmadillo

all: test_tridiag.out test_max_offdiag_symmetric.out test_jacobi_eigensolver.out\
estimate_rotations.out solve_diff.out

test_tridiag.out: tridiag.o test_tridiag.o
	$(CC) $(FLAGS) -o $@ $^

tridiag.o: tridiag.cpp
	$(CC) $(FLAGS) -c $^

test_tridiag.o: test_tridiag.cpp
	$(CC) $(FLAGS) -c $^


test_max_offdiag_symmetric.out: test_max_offdiag_symmetric.o max_offdiag_symmetric.cpp
	$(CC) $(FLAGS) -o $@ $^

test_max_offdiag_symmetric.o: test_max_offdiag_symmetric.cpp
	$(CC) $(FLAGS) -c $^

max_offdiag_symmetric.o: max_offdiag_symmetric.cpp
	$(CC) $(FLAGS) -c $^


test_jacobi_eigensolver.out: test_jacobi_eigensolver.o tridiag.o max_offdiag_symmetric.o jacobi_eigensolver.o
	$(CC) $(FLAGS) -o $@ $^

test_jacobi_eigensolver.o: test_jacobi_eigensolver.cpp
	$(CC) $(FLAGS) -c $^

jacobi_eigensolver.o: jacobi_eigensolver.cpp
	$(CC) $(FLAGS) -c $^


estimate_rotations.out: estimate_rotations.o tridiag.o jacobi_eigensolver.o max_offdiag_symmetric.o
	$(CC) $(FLAGS) -o $@ $^

estimate_rotations.o: estimate_rotations.cpp
	$(CC) $(FLAGS) -c $^


solve_diff.out: solve_diff.o tridiag.o jacobi_eigensolver.o max_offdiag_symmetric.o
	$(CC) $(FLAGS) -o $@ $^

solve_diff.o: solve_diff.cpp
	$(CC) $(FLAGS) -c $^