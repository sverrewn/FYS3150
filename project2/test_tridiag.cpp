#include <iostream>
#include <armadillo>
#include <cmath>

#include "tridiag.hpp"

void test_tridiag()
{
    int N = 6;

    arma::mat A = tridiag(N);
    
    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, A);

    check_anal(N, A, eigval, eigvec);
}

int main(int argc, char *argv[])
{
    test_tridiag();

    return 0;
    //Compile: g++ -std=c++11 test_tridiag.cpp -o test_tridiag.exe -larmadillo
    // g++ -std=c++11  test_tridiag.cpp tridiag.cpp max_offdiag_symmetric.cpp jacobi_eigensolver.cpp  -o test_tridiag.exe -larmadillo
    //Run: ./test_tridiag.exe
}