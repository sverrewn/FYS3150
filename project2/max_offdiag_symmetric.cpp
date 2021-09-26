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
                k = i;
                l = j;
            }
        }
    }

    // Return maxval 
    return maxval;
}