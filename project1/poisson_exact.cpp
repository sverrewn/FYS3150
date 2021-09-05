#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

double u(double);

int main()
{
    std::vector<double> x_vals, result;

    x_vals.reserve(1002);
    result.reserve(1002);

    for ( int i = 0; i < 1002; ++i ) {
        x_vals.push_back(i / 1001.0);
        result.push_back(u(i));
    }

    std::ofstream file;
    file.open("poisson_exact.dat");

    for ( int i = 0; i < 1002 - 1; i++ ) {
        file << x_vals[i] << ",";
    }
    file << x_vals[1001] << std::endl;

    for ( int i; i < 1002 - 1; ++i ) {
        file << std::scientific << std::setprecision(4) << result[i] << ",";
    }
    file << result[1001] << std::endl;

    return 0;
}


double u(double x)
{
    static const double u = exp(-10);
    return ( 1 - (1 - u) * x - exp(-10 * x) );
}