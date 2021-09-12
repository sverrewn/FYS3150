#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

double f(double);

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        std::cout << "Wrong amount of arguments\n" << "Usage: ./general_tridiag.out n" << std::endl;
        exit(1);
    }

    const int n = std::atoi(argv[1]);
    std::vector<double> x, v, g, temp;

    x.reserve(n); v.reserve(n); g.reserve(n), temp.reserve(n);

    const double step_size = 1.0 / (n - 1);
    for ( int i = 0; i < n; ++i ) {
        x[i] = i / ( n - 1.0 );
        g[i] = f(x[i]) * step_size * step_size;
    }

    double btemp = b[0];
    v[0] = g[0] / btemp;
    for ( int i = 1; i < n - 1; ++i ) {
        temp[i] = -1 / btemp;
        btemp = 2 + temp[i];
        v[i] = ( g[i] + v[i-1] ) / btemp;
    }

    for ( int i = n - 2; i >= 0; --i ) {
        v[i] -= temp[i+1] * v[i + 1];
    }

    std::string fname = "general_tridiag_";
    fname.append(argv[1]);
    fname.append(".dat");
    
    std::ofstream file;
    file.open(fname);

    for ( int i = 0; i < n - 1; ++i ) {
        file << x[i] << ",";
    }
    file <<x[n-1] << std::endl;

    for ( int i = 0; i < n - 1; ++i ) {
        file << std::scientific << std::setprecision(4) << v[i] << ",";
    }
    file << v[n-1] << std::endl;

    return 0;
}

double f(double x)
{
    return 100 * exp(-10*x);
}
