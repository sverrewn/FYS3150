#include <iostream>
#include <armadillo>
#include <cmath>
#include <assert.h>

#include "max_offdiag_symmetric.hpp"
#include "tridiag.hpp"
#include "jacobi_eigensolver.hpp"

int test_jacobi()
{
    int N = 6;

    // Generate random N*N matrix
    arma::mat A = tridiag(N);

    // setup
    double eps = 1e-7;

    arma::vec eigenvalues = arma::vec(N);
    arma::mat eigenvectors = arma::mat(N, N);
    
    int maxiter = 100;
    int iterations = 0;

    bool converged = false;

    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

    //test w/anal solution
    check_anal(N, A, eigenvalues, eigenvectors);

    return 0;
}

int main(int argc, char *argv[])
{
    test_jacobi();
    return 0;
    //Compile: g++ -std=c++11 test_jacobi_eigensolver.cpp -o test_jacobi_eigensolver.exe -larmadillo
    // g++ -std=c++11  test_jacobi_eigensolver.cpp tridiag.cpp max_offdiag_symmetric.cpp jacobi_eigensolver.cpp  -o test_jacobi_eigensolver.exe -larmadillo

    //Run: ./test_jacobi_eigensolver.exe
}
