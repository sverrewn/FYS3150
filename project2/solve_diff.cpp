#include <iostream>
#include <armadillo>
#include <cmath>
#include <assert.h>
#include <sstream>

#include "max_offdiag_symmetric.hpp"
#include "tridiag.hpp"
#include "jacobi_eigensolver.hpp"

void write_results(int N, arma::vec eigenvalues, arma::mat eigenvectors) 
{
    std::string fname = "data/discretization_n";
    fname.append(std::to_string(N));
    fname.append(".csv");

    std::ofstream myFile(fname);

    myFile << "# eigenvalues, eigenvectors" << std::endl; // header

    for (int i = 0; i < 3; i++) {
        myFile << eigenvalues(i) << ", ";

        for (int j = 0; j < N; j++)
        {
            if (j == N-1)
            {
                myFile << eigenvectors(i, j) << std::endl; //last line
            }

            else
            {
                myFile << eigenvectors(i, j) << ", ";
            }
        }
    }

    myFile.close();
}

void run_jacobi(int N)
{
    // Generate random tridiag N*N matrix
    arma::mat A = tridiag(N);

    double eps = 1e-7;

    arma::vec eigenvalues = arma::vec(N);
    arma::mat eigenvectors = arma::mat(N, N);

    int maxiter = 10000;
    int iterations = 0;

    bool converged = false;

    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

    arma::uvec indices = sort_index(abs(eigenvalues), "ascend");

    arma::vec smallest_eigenvalues = arma::vec(3);
    arma::mat smallest_eigenvectors = arma::mat(3, N);

    for (int i = 0; i < 3; i++)
    {
        smallest_eigenvalues(i) = eigenvalues(indices(i));
        for (int j = 0; j < N; j++)
        {
            smallest_eigenvectors(i, j) = eigenvectors(indices(j)); //vectors live in col
        }
        //smallest_eigenvectors(i) = eigenvectors(indices(i));
    }

    write_results(N, smallest_eigenvalues, smallest_eigenvectors);
}


int main(int argc, char *argv[])
{
    run_jacobi(10);
    run_jacobi(100);
    return 0;
    //Compile: g++ -std=c++11  solve_diff.cpp tridiag.cpp max_offdiag_symmetric.cpp jacobi_eigensolver.cpp  -o solve_diff.exe -larmadillo
    //Run: ./solve_diff.exe
}