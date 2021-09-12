#include <iostream>
#include <armadillo>

double lu_decomp(double);

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        std::cout << "Wrong amount of arguments\n" << "Usage: ./general_tridiag.out n" << std::endl;
        exit(1);
    }

    const int n = std::atoi(argv[1]);
    arma::vec a(n, -1), b(n, 2), c(n, -1);

    return 0;
}

int lu_decomp(arma::mat A) // A is from part 9
{
    arma::mat L, U;
    return arma::lu(L, U, A);
    
}

// Note for Mac users: Because xcode sucks - need to be run like this:  g++ -std=c++11 lu_decomposition.cpp -o lu_decomposition -larmadillo
