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

    // get the egival and eigvec, the only things we really care about
    arma::eig_sym(eigval, eigvec, A);

    std::cout << "Testing tridiag...";
    
    int errors = check_anal(N, eigval, eigvec);

    if ( errors != 0 ) {
        std::cout << "failed! Arma and analytical solution does not agree" << std::endl;
    }
    else {
        std::cout << "success! Arma and analytical solution agrees" << std::endl;
    }

    return;
}

int main()
{
    test_tridiag();

    return 0;
}