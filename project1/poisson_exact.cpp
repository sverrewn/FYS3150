#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "algorithms.hpp"

double u(double);

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        std::cout << "Wrong amount of arguments\n" << "Usage: ./poisson_exact.out n" << std::endl;
        exit(1);
    }

    const int n = std::atoi(argv[1]);

    std::vector<double> x_vals, result;

    x_vals.reserve(n);
    result.reserve(n);
    for ( int i = 0; i < n; ++i ) {
        x_vals.push_back(i / static_cast<float>(n - 1));
    }

    poisson_exact(n, x_vals, result);

    std::string fname = "poisson_exact_";
    fname.append(argv[1]);
    fname.append(".dat");

    std::ofstream file;
    file.open(fname);

    // You get horizontal columns. This makes my life easier when plotting
    for ( int i = 0; i < n - 1; i++ ) {
        file << std::scientific << std::setprecision(4) << x_vals[i] << ",";
    }
    file << x_vals[n-1] << std::endl;

    for ( int i = 0; i < n - 1; ++i ) {
        file << std::scientific << std::setprecision(4) << result[i] << ",";
    }
    file << result[n-1] << std::endl;

    return 0;
}
