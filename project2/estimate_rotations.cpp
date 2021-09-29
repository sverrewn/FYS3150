#include <iostream>
#include <armadillo>
#include <cmath>
#include <assert.h>

#include "tridiag.hpp"
#include "jacobi_eigensolver.hpp"

void estimate_rotations(const int N, int& iterations, int& maxiter, bool& converged) 
{
    // Generate random tridiag N*N matrix
    arma::mat A = tridiag(N);

    double eps = 1e-7;

    arma::vec eigenvalues = arma::vec(N);
    arma::mat eigenvectors = arma::mat(N, N);

    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
}

void write_results(int K, arma::vec n_vec, arma::vec iterations_vec, arma::vec converged_vec) 
{
    std::ofstream myFile("data/rotation_results.csv");

    myFile << "N, iterations, converged" << std::endl; // header

    for (int i = 0; i < K; i++) {
        myFile << n_vec(i) << ", " << iterations_vec(i) << ", " << converged_vec(i) << std::endl;
    }

    myFile.close();
}

void estimate_rotations_multi(int start, int stop)
{
    int K = stop - start;
    arma::vec n_vec = arma::vec(K);
    arma::vec iterations_vec = arma::vec(K);
    arma::vec converged_vec = arma::vec(K);

    int maxiter = 100000;

    int i = 0;
    for (int n = start; n < stop; n++)
    {
        int iterations = 0;
        bool converged = false;

        // get ratations
        estimate_rotations(n, iterations, maxiter, converged);

        // save values
        n_vec(i) = n;
        iterations_vec(i) = iterations;
        if (converged)
        {
            converged_vec(i) = 1;
        }

        else
        {
            converged_vec(i) = 0;
        }
        
        i++;
    }

    write_results(K, n_vec, iterations_vec, converged_vec);
}


int main()
{
    estimate_rotations_multi(6, 50);
    return 0;
}
