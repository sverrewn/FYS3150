#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

double u(double);

int main()
{
    std::vector<double> x_vals, result;

    const int n = 1002;

    x_vals.reserve(n);
    result.reserve(n);

    for ( int i = 0; i < n; ++i ) {
        x_vals.push_back(i / static_cast<float>(n - 1));
        result.push_back(u(x_vals[i]));
    }

    std::ofstream file;
    file.open("poisson_exact.dat");

    for ( int i = 0; i < n - 1; i++ ) {
        file << x_vals[i] << ",";
    }
    file << x_vals[n-1] << std::endl;

    for ( int i; i < n - 1; ++i ) {
        file << std::scientific << std::setprecision(4) << result[i] << ",";
    }
    file << result[n-1] << std::endl;

    return 0;
}


double u(double x)
{
    static const double u = exp(-10);
    return ( 1 - (1 - u) * x - exp(-10 * x) );
}