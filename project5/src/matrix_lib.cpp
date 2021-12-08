#include <armadillo>
#include <cmath>
#include <complex>


inline int translate_index(int i, int j, int length)
{
    return ( i * length + j);
}


void create_a(arma::cx_vec& a, int n, int dt, double r, arma::cx_mat& v)
{
    for ( int i = 0; i < n; ++i ) {
        for ( int j = 0; j < n; ++j ) {
            a.at(translate_index(i,j, n)) = 1 + 4*r + ( arma::cx_double(0,dt / 2 ) * v.at(i,j));
        }
    }

    return;
}


void create_b(arma::cx_vec& b, int n, int dt, double r, arma::cx_mat& v)
{
    for ( int i = 0; i < n; ++i ) {
        for ( int j = 0; j < n; ++j )
            b.at(translate_index(i,j, n)) = 1 - 4*r + ( arma::cx_double(0, dt / 2 ) * v.at(i,j));
    }

    return;
}


void initial_u(arma::cx_vec& u, int len_x, int len_y, double x_c, double y_c, double p_x, double p_y, double sig_x, double sig_y)
{   
    for ( int k = 0; k < len_x; ++k ) {
        for ( int l = 0; l < len_y; ++l ) {
            u.at(translate_index(k,l, len_x)) = std::exp( 
                -( (k-x_c) * (k-x_c) ) / ( 2 * sig_x * sig_x )
                -( (l-y_c) * (l-y_c) ) / ( 2 * sig_y * sig_y )
                + arma::cx_double(0, p_x * (k-x_c)) + arma::cx_double(0, p_y * (l-y_c))
            );
        }
    }


    arma::cx_vec u_star = arma::conj(u);

    arma::cx_double sum = 0.;    

    for ( int i = 0; i < u.size(); ++i ) {
        sum += u_star[i] * u[i];
    }

    u = u / std::sqrt(sum);

    return;
}


void init_potential(arma::mat V, int M)
{
    double x_w = 0.02, wall_pos = 0.5;
    double slit_dist = 0.005, slit_op = 0.05;

    int wall_width = std::round(x_w * M);
    int wall_center = std::round(wall_pos*M);
    int slit_distance = std::round(slit_dist * M);
    int slit_opening = std::round(slit_op * M);

    for ( int i = 0; i < M; ++i ) {
        V.at(0, i) = 10000;
        V.at(i, 0) = 10000;
        V.at(i, M-1) = 10000;
        V.at(M-1, i) = 10000;
    }

    for ( int j = 0; j < M; ++j ) {

    }

    int wall_start = wall_center - std::round(wall_width / 2);

    int slit1_start = wall_center - std::round(slit_distance/2) - slit_opening;
    int slit1_end = wall_center - std::round(slit_distance/2);

    int slit2_start = wall_center + std::round(slit_distance/2);
    int slit2_end = wall_center - std::round(slit_distance/2) + slit_opening;

    for ( int i = 0; i < M; ++i ) { // rows
        for ( int k = wall_start; k <= wall_start + wall_width; ++k ) {
            V.at(i,k) = 1;
        }
    }

    for ( int i = slit1_start; i <= slit1_end; ++i ) {
        for ( int k = wall_start; k <= wall_start + wall_width; ++k ) {
            V.at(i, k) = 1;
        }
    }

    for ( int i = slit2_start; i <= slit2_end; ++i ) {
        for ( int k = wall_start; k <= wall_start + wall_width; ++k ) {
            V.at(i, k) = 1;
        }
    }

    return;
}


void fill_matrices(arma::cx_mat& A, arma::cx_mat& B, double r, arma::cx_vec& a, arma::cx_vec& b, int sub_len)
{
    int tot_len = sub_len * sub_len;
    
    /* ### CREATE A ### */

    // top row matrices
    // middle matrix
    A.at(0,0) = a[0]; 
    A.at(0,1) = -r;
    
    // right matrix
    A.at(0, sub_len) = -r;

    for ( int i = 1; i < sub_len - 1; ++i ) {
        // middle matrix
        A.at(i, i-1) = -r;
        A.at(i, i)   = a[i];
        A.at(i, i+1) = -r;
        
        // right matrix
        A.at(i, sub_len + i) = -r;
    }
    A.at(sub_len - 1, sub_len-2) = -r;
    A.at(sub_len - 1, sub_len - 1) = a[sub_len-1];
    A.at(sub_len - 1, 2*sub_len - 1) = -r;

    // all middle-row matrices
    for ( int i = 1; i < sub_len - 1; ++i ) {
        // "(0,0)" for matrix on diagonal. offset by sub_len to get the other matrices of a given row
        int idx = i * sub_len;
        
        // Top row of sub matrices
        // left matrix
        A.at(idx, idx - sub_len) = -r;
        
        // middle matrix
        A.at(idx, idx) = a[idx];
        A.at(idx, idx + 1) = -r;

        // right matrix
        A.at(idx, idx + sub_len) = -r;

        // middle parts of the matrices
        for ( int k = idx + 1; k < idx + sub_len - 1; ++k ) {
            // left matrix
            A.at(k, k - sub_len) = -r;

            // middle matrix
            A.at(k, k-1) = -r;
            A.at(k, k) = a[k];
            A.at(k, k+1) = -r;

            // right matrix
            A.at(k, k + sub_len) = -r;
        }
        // Last row of sub matrices
        // left matrix
        int l_idx = idx + sub_len - 1;
        A.at(l_idx, l_idx - sub_len) = -r;
        
        // middle matrix
        A.at(l_idx, l_idx - 1) = -r;
        A.at(l_idx, l_idx) = a[l_idx];

        // right matrix
        A.at(l_idx, l_idx + sub_len) = -r;
    }
    // bottom row
    int idx = (sub_len - 1) * sub_len;
    // left matrix
    A.at(idx, idx - sub_len) = -r;

    // middle matrix
    A.at(idx, idx) = a[idx];
    A.at(idx, idx + 1) = -r;

    for ( int i = idx + 1; i < tot_len - 1; ++i ) {
        // left matrix
        A.at(i, i - sub_len) = -r;
        
        // middle matrix
        A.at(i, i-1) = -r;
        A.at(i, i) = a[i];
        A.at(i, i+1) = -r;
    }
    // left matrix
    A.at(tot_len - 1, tot_len-1 - sub_len) = -r;

    // middle matrix
    A.at(tot_len - 1, tot_len - 2) = -r;
    A.at(tot_len - 1, tot_len - 1) = a[tot_len - 1];

    /* ### CREATE B ### */

    // top row matrices
    // middle matrix
    B.at(0,0) = b[0]; 
    B.at(0,1) = r;
    
    // right matrix
    B.at(0, sub_len) = r;

    for ( int i = 1; i < sub_len - 1; ++i ) {
        // middle matrix
        B.at(i, i-1) = r;
        B.at(i, i)   = b[i];
        B.at(i, i+1) = r;
        
        // right matrix
        B.at(i, sub_len + i) = r;
    }
    B.at(sub_len - 1, sub_len-2) = r;
    B.at(sub_len - 1, sub_len - 1) = b[sub_len-1];
    B.at(sub_len - 1, 2*sub_len - 1) = r;

    // all middle-row matrices
    for ( int i = 1; i < sub_len - 1; ++i ) {
        // "(0,0)" for matrix on diagonal. offset by sub_len to get the other matrices of a given row
        int idx = i * sub_len;
        
        // Top row of sub matrices
        // left matrix
        B.at(idx, idx - sub_len) = r;
        
        // middle matrix
        B.at(idx, idx) = b[idx];
        B.at(idx, idx + 1) = r;

        // right matrix
        B.at(idx, idx + sub_len) = r;

        // middle parts of the matrices
        for ( int k = idx + 1; k < idx + sub_len - 1; ++k ) {
            // left matrix
            B.at(k, k - sub_len) = r;

            // middle matrix
            B.at(k, k-1) = r;
            B.at(k, k) = b[k];
            B.at(k, k+1) = r;

            // right matrix
            B.at(k, k + sub_len) = r;
        }
        // Last row of sub matrices
        // left matrix
        int l_idx = idx + sub_len - 1;
        B.at(l_idx, l_idx - sub_len) = r;
        
        // middle matrix
        B.at(l_idx, l_idx - 1) = r;
        B.at(l_idx, l_idx) = b[l_idx];

        // right matrix
        B.at(l_idx, l_idx + sub_len) = r;
    }
    // bottom row
    idx = (sub_len - 1) * sub_len;
    // left matrix
    B.at(idx, idx - sub_len) = r;

    // middle matrix
    B.at(idx, idx) = b[idx];
    B.at(idx, idx + 1) = r;

    for ( int i = idx + 1; i < tot_len - 1; ++i ) {
        // left matrix
        B.at(i, i - sub_len) = r;

        // middle matrix
        B.at(i, i-1) = r;
        B.at(i, i) = b[i];
        B.at(i, i+1) = r;
    }
    // left matrix
    B.at(tot_len - 1, tot_len-1 - sub_len) = r;

    // middle matrix
    B.at(tot_len - 1, tot_len - 2) = r;
    B.at(tot_len - 1, tot_len - 1) = b[tot_len - 1];

    return;
}


void solve_eqs(arma::cx_mat A, arma::cx_mat& B, arma::cx_vec b, arma::cx_vec& u)
{
    b = B * u;
    u = arma::solve(A, b);

    return;
}