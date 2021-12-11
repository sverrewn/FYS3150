#include <armadillo>
#include <cmath>
#include <complex>


// convert matrix index i,j to index u for a vector consisting of all columns of said matrix
inline int translate_index(int i, int j, int length)
{
    return ( i * length + j);
}


// Requires a to be preallocated with a size of (M-2)^2
// Requires v (the potential) to be initialised
void init_a(arma::cx_vec& a, int n, int dt, arma::cx_double r, arma::mat& v)
{
    for ( int i = 0; i < n; ++i ) {
        for ( int j = 0; j < n; ++j ) {
            a.at(translate_index(i,j, n)) = 1. + 4. * r + ( arma::cx_double(0,dt / 2 ) * v.at(i,j));
        }
    }

    return;
}


// Requires b to be preallocated with a size of (M-2)^2
// Requires v (the potential) to be initialised
void init_b(arma::cx_vec& b, int n, int dt, arma::cx_double r, arma::mat& v)
{
    for ( int i = 0; i < n; ++i ) {
        for ( int j = 0; j < n; ++j )
            b.at(translate_index(i,j, n)) = 1. - 4. * r + ( arma::cx_double(0, dt / 2 ) * v.at(i,j));
    }

    return;
}


// Sets up u for t=0. Requires u to be preallocated with a size of (M-2)^2
void initial_u(arma::cx_vec& u, int len_x, int len_y, double x_c, double y_c, double p_x, double p_y, double sig_x, double sig_y, double h)
{   
    // y value
    for ( int l = 1; l <= len_y; ++l ) {
        // x value
        for ( int k = 1; k <= len_x; ++k ) {
            u.at(translate_index(k-1,l-1, len_x)) = std::exp( 
                -( (k*h - x_c) * (k*h - x_c) ) / ( 2 * sig_x * sig_x )
                -( (l*h - y_c) * (l*h - y_c) ) / ( 2 * sig_y * sig_y )
                + arma::cx_double(0, p_x * (k*h-x_c)) + arma::cx_double(0, p_y * (l*h-y_c))
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


// Sets up the wall and the slits. Requires V to be preallocated with a size of M and set to 0
void init_potential(arma::mat& V, double v, int M, int slits)
{   
    if ( slits < 1 ) { // no wall needed
        std::cout << "NO SLITS" << std::endl;
        return;
    }

    double x_w = 0.02, wall_pos = 0.5;
    double slit_dist = 0.05, slit_op = 0.05;

    // x-plane values for the wall
    int wall_width = std::round(x_w * M);
    int wall_center = std::round(wall_pos * M);
    int wall_start = wall_center - std::round(wall_width / 2);
    
    int slit_distance = std::round(slit_dist * M); // distance between slits
    int slit_opening = std::round(slit_op * M); // slit length

    // Set up wall in the middle
    for ( int i = 0; i < M; ++i ) { // rows
        for ( int k = wall_start; k <= wall_start + wall_width; ++k ) {
            V.at(i,k) = v;
            std::cout << i << ", " << k << std::endl;
        }
    }
    
   // std::cout << V << std::endl;

    int cent_y = std::round(M/2);

    if ( slits % 2 == 1 ) { // odd number of slits, meaning one is in the center
        std::cout << "Creating " << slits << " slits" << std::endl;
        // there's a slit in the middle
        int center_offset = std::round(slit_opening / 2);
        // making the slits going from the center and outwards in both directions simultaneously
        for ( int i = 1; i <= (( slits - 1) / 2 ); ++i ) {
            for ( int k = wall_start; k <= wall_start + wall_width; ++k ) {
                for ( int j = 0; j < slit_distance; ++j ) {
                    // offset taking into account there's a slit in the middle
                    int offset = center_offset + i * slit_distance + (i-1) * slit_opening;
                    V.at(cent_y + offset + j, k) = 0;
                    V.at(cent_y - offset - j, k) = 0;
                }
            }
        }

        // drawing the center slit
        for ( int i = cent_y - center_offset; i <= cent_y + center_offset; ++i ) {
            for ( int k = wall_start; k <= wall_start + wall_width; ++k ) {
                V.at(i, k) = 0;
            }
        }
    }
    else { // even number slits, no opening in the middle
        std::cout << "Creating " << slits << " slits" << std::endl;
        // no slit in the middle
        int center_offset = std::round(slit_distance / 2);

        // making the slits going from the center and outwards in both directions simultaneously
        for ( int i = 0; i < ( slits / 2); ++i ) {
            for ( int k = wall_start; k <= wall_start + wall_width; ++k ) {
                for ( int j = 0; j < slit_distance; ++j ) {
                    // no slit in the middle, slits begin after center_offset
                    int offset = center_offset + i * (slit_distance + slit_opening);
                    V.at(cent_y + offset + j, k) = 0;
                    V.at(cent_y - offset - j, k) = 0;
                    std::cout << cent_y + offset + j << ", " << k << std::endl;
                }
            }
        }

        std::cout << V << std::endl;
    }

    return;
}


// Fills the matrix A and B. 
// Requires A and B to be preallocated with a size of (M-2)^2 x (M-2)^2
// Has to be called after a and b are initalised
void fill_matrices(arma::sp_cx_mat& A, arma::sp_cx_mat& B, arma::cx_double r, arma::cx_vec& a, arma::cx_vec& b, int sub_len)
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


// Solves the equation Au^(n+1) = Bu^n for u^(n+1)
void solve_eqs(arma::sp_cx_mat& A, arma::sp_cx_mat& B, arma::cx_vec& b, arma::cx_vec& u)
{
    b = B * u;
    u = arma::spsolve(A, b);

    return;
}