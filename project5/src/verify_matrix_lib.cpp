#include <armadillo>
#include <iostream>
#include <vector>

#include "matrix_lib.hpp"


void test_matrix();
void test_potential();


int main()
{   
    test_matrix();
    //test_potential();

    return 0;
}


void test_matrix()
{
    int size = 9;
    int sub_size = 3;

    arma::sp_cx_mat A = arma::sp_cx_mat(size, size);
    arma::sp_cx_mat B = arma::sp_cx_mat(size, size);

    double r = 1;

    arma::cx_vec a = {0,10,20,30,40,50,60,70,80};
    arma::cx_vec b = {80,70,60,50,40,30,20,10,0};

    fill_matrices(A, B, r, a, b, sub_size);

    arma::cx_mat C = arma::cx_mat(A);
    arma::cx_mat D = arma::cx_mat(B);
    std::cout << C << std::endl;
    std::cout << "\n" << D << std::endl;

    size = 16; sub_size = 4;

    A = arma::sp_cx_mat(size, size);
    B = arma::sp_cx_mat(size, size);

    a = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    b = {15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0};

    fill_matrices(A,B,r,a,b,sub_size);

    C = arma::cx_mat(A);
    D = arma::cx_mat(B);
    std::cout << C << std::endl;
    std::cout << D << std::endl;
}


void test_potential()
{   
    double h = 0.01;
    int len = 1 / h + 1;
    double v_0 = 5.;
    int slits = 2;

    arma::mat V = arma::mat(len, len, arma::fill::zeros);
    init_potential(V, v_0, len, slits);

    std::cout << V << std::endl;
}