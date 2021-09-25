#include <iostream>
#include <armadillo>
#include <cmath>
#include <assert.h>

// Determine the the max off-diagonal element of a symmetric matrix A
// - Saves the matrix element indicies to k and l 
// - Returns absolute value of A(k,l) as the function return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l)
{
    // Get size of the matrix A
    int N = A.n_rows;  

    // Initialize references k and l to the first off-diagonal element of A
    k = 0;
    l = 1;

    // Initialize a double variable 'maxval' to A(k,l). We'll use this variable to keep track of the largest off-diag element.
    double maxval = A(k, l);

    // Loop through all elements in the upper triangle of A (not including the diagonal)
    // When encountering a matrix element with larger absolute value than the current value of maxval,
    // update k, l and max accordingly.
    for (int i = 0; i < N-1; i++) 
    {
        for (int j = i + 1; j < N; j++) 
        {
            if (abs(A(i,j)) > maxval)
            {
                maxval = abs(A(i,j));
            }
        }
    }

    // Return maxval 
    return maxval;
}

int test_max_offdiag_symmetric() {

    // creates test matrix
    int N = 4;

    arma::mat A(N, N, arma::fill::eye); // eye: off-diag = 0, diag = 1
    A(0, 3) = A(3, 0) = 0.5;
    A(1, 3) = A(3, 1) = -0.7;

    //const arma::mat B = A;

    //creates k & l
    int k, l;

    //tests max_offdiag_symmetric()
    float maxval = max_offdiag_symmetric(A, k, l);

    // assert test
    float eps = 1e-7;
    assert(abs(abs(maxval) - 0.7) < eps);

    return 0;
}

int main(int argc, char *argv[])
{
    test_max_offdiag_symmetric(); //call testfunction
    return 0;
    //Compile: g++ -std=c++11 max_offdiag_symmetric.cpp -o max_offdiag_symmetric.exe -larmadillo
    //Run: ./max_offdiag_symmetric.exe
}