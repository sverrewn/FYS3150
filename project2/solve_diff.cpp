#include <iostream>
#include <armadillo>
#include <cmath>
#include <assert.h>
#include <sstream>

#include "max_offdiag_symmetric.hpp"
#include "tridiag.hpp"
#include "jacobi_eigensolver.hpp"

void write_results(int N, arma::vec eigenvalues, arma::mat eigenvectors, bool anal=false) 
{
    std::string fname;

    if (anal)
    {
        fname = "data/anal_n";
    }

    else
    {
        fname = "data/discretization_n";
    }
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

void run_anal(int N, arma::uvec indices)
{
    int a = -1;
    int d = 2;

    // creates eigval w/analytical formula
    arma::vec anal_eigval(N);

    for (int i = 1; i <= N; i++)
    {
        anal_eigval[i-1] = d + 2 * a * cos((i * M_PI) / (N + 1));
    }

    // creates eigvec w/analytical formula
    arma::mat anal_eigvec(N, N);

    for (int i = 1; i <= N; i++) 
    {
        for (int n = 1; n <= N; n++) 
        {
            anal_eigvec.at(i-1, n-1) = (sin((n * i * M_PI) / (N + 1)));
        }
    }

    anal_eigvec.t();
    
    arma::mat anal_eigvec_normalise = normalise(anal_eigvec);

    arma::vec anal_smallest_eigenvalues = arma::vec(3);
    arma::mat anal_smallest_eigenvectors = arma::mat(3,N);

    for (int i = 0; i < 3; i++)
    {
        anal_smallest_eigenvalues(i) = anal_eigval(indices(i));
        for (int j = 0; j < N; j++)
        {
            anal_smallest_eigenvectors(i,j) = anal_eigvec_normalise(j,indices(i)); //vectors live in col
        }
    }

    write_results(N, anal_smallest_eigenvalues, anal_smallest_eigenvectors, true);
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

    eigenvalues.print("eigenvalues=");
    eigenvectors.print("eigenvectors=");

    arma::uvec indices = sort_index(abs(eigenvalues), "ascend");

    arma::vec smallest_eigenvalues = arma::vec(3);
    arma::mat smallest_eigenvectors = arma::mat(3,N);

    for (int i = 0; i < 3; i++)
    {
        smallest_eigenvalues(i) = eigenvalues(indices(i));
        for (int j = 0; j < N; j++)
        {
            smallest_eigenvectors(i,j) = eigenvectors(j,indices(i)); //vectors live in col
        }
        //smallest_eigenvectors(i) = eigenvectors(indices(i));
    }

    smallest_eigenvectors.print("smallest_eigenvectors=");

    write_results(N, smallest_eigenvalues, smallest_eigenvectors);

    run_anal(N, indices);
}

int main(int argc, char *argv[])
{
    //N = 10
    run_jacobi(10);

    // N = 100
    run_jacobi(100);
    
    return 0;
    //Compile: g++ -std=c++11  solve_diff.cpp tridiag.cpp max_offdiag_symmetric.cpp jacobi_eigensolver.cpp  -o solve_diff.exe -larmadillo
    //Run: ./solve_diff.exe
}