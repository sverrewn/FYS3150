#include <iostream>
#include <armadillo>
#include <cmath>
#include "tridiag.hpp"

arma::mat tridiag(int N) // Creates a NxN tridiag matrix
{
    arma::mat A(N, N, arma::fill::zeros);

    A.diag(-1) += -1;
    A.diag(0) += 2;
    A.diag(1) += -1; 

    return A;
}

int check_anal(int N, arma::mat A, arma::vec eigval, arma::mat eigvec) //Checks with analytical solution
{
    double h = 1.0 / (N + 1.0);
    double a = -1.0 / pow(h, 2);
    double d = 2.0 / pow(h, 2);

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

    float eps = exp(-7); // = 0

    // check if it corresponds with analytical solution
    if ((arma::sum(abs(eigval - anal_eigval))) > eps) 
    {
        throw std::runtime_error("Eigenvalues does not correspond with analytical solution");
    }


    for (int i = 0; i < N; i++) 
    {
        if ((arma::sum(abs(eigvec[i] + anal_eigvec_normalise[i]))) > eps) 
        {
            throw std::runtime_error("Eigenvec does not correspond with analytical solution");
        }
    }

    return 0;
} 





