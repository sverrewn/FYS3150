#include <armadillo>
#include <iostream>
#include <vector>

#include "matrix_lib.hpp"

int main()
{
    int size = 9;
    int sub_size = 3;

    arma::cx_mat A = arma::cx_mat(size, size, arma::fill::zeros);
    arma::cx_mat B = arma::cx_mat(size, size, arma::fill::zeros);

    double r = 1;

    arma::cx_vec a = {0,10,20,30,40,50,60,70,80};
    arma::cx_vec b = {80,70,60,50,40,30,20,10,0};

    fill_matrices(A, B, r, a, b, sub_size);

    std::cout << A << std::endl;
    std::cout << "\n" << B << std::endl;
}