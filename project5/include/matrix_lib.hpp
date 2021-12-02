#ifndef SNIPPETS_HPP
#define SNIPPETS_HPP

int translate_index(int i, int j, int length);
void fill_matrices(arma::mat& A, arma::mat& B, double r, arma::vec& a, arma::vec& b, int sub_len);
void solve_eqs(arma::mat& A, arma::mat& B, arma::vec& b, arma::vec& u);
#endif