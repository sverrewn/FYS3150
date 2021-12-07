#ifndef SNIPPETS_HPP
#define SNIPPETS_HPP

int translate_index(int i, int j, int length);
void fill_matrices(arma::cx_mat& A, arma::cx_mat& B, double r, arma::cx_vec& a, arma::cx_vec& b, int sub_len);
void solve_eqs(arma::cx_mat& A, arma::cx_mat& B, arma::cx_vec& b, arma::cx_vec& u);
#endif