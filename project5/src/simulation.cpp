#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "matrix_lib.hpp"

void u_mat(int len, arma::cx_vec& u, arma::cx_mat& m)
{
    for (int i=0; i < len; i++){
        for (int j=0; j < len; j++){
            m.at(j,i) = u(translate_index(i, j, len));
        }
    }
    return;
}

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

    int len = 1 / h - 1;
    arma::mat V = arma::mat(len, len, arma::fill::zeros);
    init_potential(V, v_0, len, slits);
    V.save("data/potential.bin");

    int len_sq = (len) * (len);

    arma::cx_vec u = arma::cx_vec(len_sq);
    initial_u(u, len + 2, len + 2, x_c, y_c, p_x, p_y, sig_x, sig_y, h);

    arma::cx_double r = arma::cx_double(0, dt / (2 * h * h));

    arma::cx_vec a = arma::cx_vec(len_sq);
    arma::cx_vec b = arma::cx_vec(len_sq);
    init_a(a, len, dt, r, V);
    init_b(b, len, dt, r, V);


    arma::sp_cx_mat A = arma::sp_cx_mat(len_sq, len_sq);
    arma::sp_cx_mat B = arma::sp_cx_mat(len_sq, len_sq);
    fill_matrices(A, B, r, a, b, len);

    
    int len_T = T/dt + 1;
    int i = 0;
    arma::cx_cube U = arma::cx_cube(len, len, len_T);
    arma::cx_mat U_i = arma::cx_mat(len, len);
    u_mat(len, u, U_i);
    U.slice(i++) = U_i;

    std::cout << T << " " << len_T << std::endl;
    for (double t = dt; t <= T; t += dt){
        solve_eqs(A, B, b, u);
        u_mat(len, u, U_i);
        U.slice(i++) = U_i;
        std::cout << t << std::endl;
    }

    U.save(fname);
    
    return 0;
}