#include <armadillo>


inline int translate_index(int i, int j, int length)
{
    return ( i * length + j);
}


void fill_matrices(arma::mat& A, arma::mat& B, double r, arma::vec& a, arma::vec& b, int sub_len)
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
}


void solve_eqs(arma::mat A, arma::mat& B, arma::vec b, arma::vec& u)
{
    b = B * u;
    u = arma::solve(A, b);
}