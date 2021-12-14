#ifndef MATRIX_LIB_HPP
#define MATRIX_LIB_HPP

int translate_index(int i, int j, int length);
void init_a(arma::cx_vec& a, int n, double dt, arma::cx_double r, arma::mat& v);
void init_b(arma::cx_vec& b, int n, double dt, arma::cx_double r, arma::mat& v);
void initial_u(arma::cx_vec& u, int len_x, int len_y, double x_c, double y_c, double p_x, double p_y, double sig_x, double sig_y, double h);
void init_potential(arma::mat& V, double v, int M, int slits);
void fill_matrices(arma::sp_cx_mat& A, arma::sp_cx_mat& B, arma::cx_double r, arma::cx_vec& a, arma::cx_vec& b, int sub_len);
void solve_eqs(arma::sp_cx_mat& A, arma::sp_cx_mat& B, arma::cx_vec& b, arma::cx_vec& u);

#endif