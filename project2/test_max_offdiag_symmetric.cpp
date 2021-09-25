#include <iostream>
#include <armadillo>
#include <cmath>
#include <assert.h>

#include "max_offdiag_symmetric.cpp"

int test_max_offdiag_symmetric() {

    // creates test matrix
    int N = 4;

    arma::mat A(N, N, arma::fill::eye); // eye: off-diag = 0, diag = 1
    A(0, 3) = A(3, 0) = 0.5;
    A(1, 3) = A(3, 1) = -0.7;

    //const arma::mat B = A;

    //creates k & l
    int k, l;

    //tests max_offdiag_symmetric()
    float maxval = max_offdiag_symmetric(A, k, l);

    // assert test
    float eps = 1e-7;
    assert(abs(abs(maxval) - 0.7) < eps);

    return 0;
}

int main(int argc, char *argv[])
{
    test_max_offdiag_symmetric(); //call testfunction
    return 0;
    //Compile: g++ -std=c++11 max_offdiag_symmetric.cpp -o max_offdiag_symmetric.exe -larmadillo
    //Run: ./max_offdiag_symmetric.exe
}