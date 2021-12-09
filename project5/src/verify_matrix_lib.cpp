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

    size = 16; sub_size = 4;

    A = arma::cx_mat(size, size, arma::fill::zeros);
    B = arma::cx_mat(size, size, arma::fill::zeros);

    a = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    b = {15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0};

    fill_matrices(A,B,r,a,b,sub_size);

    std::cout << A << std::endl;
    std::cout << B << std::endl;

    return 0;
}