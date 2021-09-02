#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

double u(double);

int main()
{
    std::vector<double> x_vals, result;

    for ( int i = 0; i < 1002; ++i ) {
        x_vals.push_back(i / 1001.0);
    }

    for ( double x: x_vals ) {
        result.push_back(u(x));
    }

    std::ofstream file;
    file.open("poisson_exact.dat");

    for ( double x: x_vals) {
        file << x << ",";
    }
    file << std::endl;
    for ( double x: result) {
        file << std::scientific << std::setprecision(4) << x << ",";
    }
    file << std::endl;

    return 0;
}

double u(double x)
{
    return ( 1 - (1 - exp(-10)) * x - exp(-10 * x) );
}