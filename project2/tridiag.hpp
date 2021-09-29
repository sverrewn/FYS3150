#ifndef TRIDIAG_HPP
#define TRIDIAG_HPP

arma::mat tridiag(int N);

int check_anal(int N, arma::vec eigval, arma::mat eigvec);

#endif