#include <iostream>
#include <armadillo>
#include <cmath>

#include "max_offdiag_symmetric.hpp"

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

    int error {0};
    
    if (std::abs(maxval - 0.7) > eps) {
        error = 1;
        std::cout << "maxval: " << maxval << std::endl;
        std::cout << "k: " << k << " l: " << l << std::endl;
    }

    return error;
}

int main()
{   
    std::cout << "Testing max_offdiag_symmetric...";
    int error = test_max_offdiag_symmetric(); //call testfunction
    
    if (error != 0 ) {
        std::cout << "failed!" << std::endl;
    }
    else {
        std::cout << "success!" << std::endl;
    }

    return 0;
}