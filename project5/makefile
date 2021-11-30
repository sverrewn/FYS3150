CC = g++ -std=c++11
FLAGS = -Wextra -Wpedantic -Wall -I include/

test: verify_snippet.out

verify_snippet.out: verify_snippet.o snippets.o
	$(CC) $(FLAGS) -fopenmp -o $@ $^

verify_snippet.o: src/verify_snippet.cpp
	$(CC) $(FLAGS) -c $^

snippets.o: src/snippets.cpp
	$(CC) $(FLAGS) -c $^