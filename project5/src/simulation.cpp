#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "matrix_lib.hpp"


int main(int argc, char* argv[]) 
{   
    if ( argc < 13 ) {
        std::cout << "Missing command line arguments" << std::endl;

        return 1;
    }
    // h, dt, T, x_c, sig_x, p_x, y_c, sig_y, p_y, v_0
    double h = atof(argv[1]);
    double dt = atof(argv[2]);
    double T = atof(argv[3]);
    double x_c = atof(argv[4]);
    double sig_x = atof(argv[5]);
    double p_x = atof(argv[6]);
    double y_c = atof(argv[7]);
    double sig_y = atof(argv[8]);
    double p_y = atof(argv[9]);
    double v_0 = atof(argv[10]);
    int slits = atoi(argv[11]);
    std::string fname = argv[12];

    int len = 1 / h + 1;
    arma::mat V = arma::mat(len, len, arma::fill::zeros);
    init_potential(V, v_0, len, slits);

    int len_sq = (len - 2) * (len - 2);

    arma::cx_vec u = arma::cx_vec(len_sq);
    initial_u(u, len, len, x_c, y_c, p_x, p_y, sig_x, sig_y, h);

    arma::cx_double r = arma::cx_double(0, dt / (2 * h * h));

    arma::cx_vec a = arma::cx_vec(len_sq);
    arma::cx_vec b = arma::cx_vec(len_sq);
    init_a(a, len-2, dt, r, V);
    init_b(b, len-2, dt, r, V);


    arma::sp_cx_mat A = arma::sp_cx_mat(len_sq, len_sq);
    arma::sp_cx_mat B = arma::sp_cx_mat(len_sq, len_sq);
    fill_matrices(A, B, r, a, b, len - 2);

    std::vector<arma::cx_vec> u_hist;
    u_hist.push_back(u);

    for ( double t = dt; t <= T; t += dt ) {
        solve_eqs(A, B, b, u);
        u_hist.push_back(u);
        std::cout << t << std::endl;
    }

    std::ofstream file;
    file.open(fname);

    for (auto f : u_hist ) {
        for ( auto v : f ) {
            file << v << " ";
        }
        file << std::endl;
    }

    return 0;
}