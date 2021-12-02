#include <armadillo>
#include <iostream>
#include <vector>

#include "matrix_lib.hpp"

int main()
{
    int size = 9;
    int sub_size = 3;

    arma::mat A = arma::mat(size, size, arma::fill::zeros);
    arma::mat B = arma::mat(size, size, arma::fill::zeros);

    double r = 1;

    arma::vec a = {0,10,20,30,40,50,60,70,80};
    arma::vec b = {80,70,60,50,40,30,20,10,0};

    fill_matrices(A, B, r, a, b, sub_size);

    std::cout << A << std::endl;
    std::cout << "\n" << B << std::endl;
}