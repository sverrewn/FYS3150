#ifndef SNIPPETS_HPP
#define SNIPPETS_HPP

int translate_index(int i, int j, int length);
void create_a(arma::cx_vec& a, int n, int dt, double r, arma::cx_mat& v);
void create_b(arma::cx_vec& b, int n, int dt, double r, arma::cx_mat& v);
void initial_u(arma::cx_vec& u, int len_x, int len_y, double x_c, double y_c, double p_x, double p_y, double sig_x, double sig_y);
void init_potential(arma::mat V, int M);
void fill_matrices(arma::cx_mat& A, arma::cx_mat& B, double r, arma::cx_vec& a, arma::cx_vec& b, int sub_len);
void solve_eqs(arma::cx_mat& A, arma::cx_mat& B, arma::cx_vec& b, arma::cx_vec& u);
#endif