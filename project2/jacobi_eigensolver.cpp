#include <iostream>
#include <armadillo>
#include <cmath>
#include <assert.h>

#include "max_offdiag_symmetric.hpp"
#include "tridiag.hpp"

using namespace std;
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l) 
{
    // Performs a single Jacobi rotation, to "rotate away"
    // the off-diagonal element at A(k,l).
    // - Assumes symmetric matrix, so we only consider k < l
    // - Modifies the input matrices A and R

    // Code from textbook p. 219-220

    float s, c; //sin, cos
    int N = A.n_cols;

    // Performs a single Jacobi rotation, to "rotate away"
    if ((A(k, l) != 0.0)) 
    {
        float t, tau;
        
        tau = (A(l, l) - A(k, k))/(2*A(k, l));

        if (tau > 0) 
        {
            t = 1.0 / (tau + sqrt(1.0 + tau * tau));
        }

        else 
        {
            t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
        }

        c = 1 / sqrt(1+t*t);
        s = c * t;
    }

    else
    {
        c = 1.0;
        s = 0.0;
    }
    
    float A_kk, A_ll, A_ik, A_il, R_ik, R_il;

    A_kk = A(k, k);
    A_ll = A(l, l);

    //change all elements indexed at somethin with k or l
    A(k, k) = c * c * A_kk - 2.0 * c * s * A(k, l) + s * s * A_ll;
    A(l, l) = s * s * A_kk + 2.0 * c * s * A(k, l) + c * c * A_ll;
    A(k,l) = A(l, k) = 0.0;

    //changing remaining elements
    for (int i = 0; i < N; i++) 
    {
        if (i != k && i != l)
        {
            A_ik = A(i, k);
            A_il = A(i, l);

            A(i, k) = c * A_ik - s * A_il;
            A(k, i) = A(i, k);

            A(i, l) = c * A_il + s * A_ik;
            A(l, i) = A(i, l);
        }

        //compute new eigvecs
        R_ik = R(i, k);
        R_il = R(i, l);

        R(i, k) = c * R_ik - s * R_il;
        R(i, l) = c * R_il + s * R_ik;
    }
}


void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged)
{
    // Jacobi method eigensolver:
    // - Runs jacobo_rotate until max off-diagonal element < eps
    // - Writes the eigenvalues as entries in the vector "eigenvalues"
    // - Writes the eigenvectors as columns in the matrix "eigenvectors"
    //   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
    // - Stops if it the number of iterations reaches "maxiter"
    // - Writes the number of iterations to the integer "iterations"
    // - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter

    arma::mat R = arma::mat(arma::size(A));
    arma::mat A_copy = A;

    int k = 0, l = 0;

    float max_off_diag = 10;

    while ((max_off_diag > eps) && (iterations < maxiter))
    {
        max_off_diag = max_offdiag_symmetric(A_copy, k, l);
        try
        {
            jacobi_rotate(A_copy, R, k, l);
        }
        catch(const std::exception& e)
        {
            std::cout << "failed jacobi" << max_off_diag  << " "  <<eps  << '\n';
        }
        try
        {
            arma::eig_sym(eigenvalues, eigenvectors, A);
        }
        catch(const std::exception& e)
        {
            std::cout << "failed eig_sym" << max_off_diag  << " "  <<eps  << '\n';
        }

        iterations++;
    }

    if (iterations < maxiter) 
    {
        converged = true;
        //std::cout << "Converged!" << std::endl;
    }
} 

