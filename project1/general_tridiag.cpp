#include <vector>
#include <cmath>
#include <iostream>

double f(double);

int main(int argc, char *argv[]) {
    if ( argc != 2 ) {
        std::cout << "Wrong amount of arguments\n" << "Usage: ./general_tridiag n" << std::endl;
    }
    const int n = std::atoi(argv[1]);
    std::vector<double> a(n, -1), b(n, 2), c(n, -1), x, v, g, bt, gt, temp;

    x.reserve(n); v.reserve(n); g.reserve(n), bt.reserve(n), gt.reserve(n), temp.reserve(n);

    const double step_size = 1 / (n - 1.0);
    for ( int i = 0; i < n; ++i ) {
        x[i] = i / ( n - 1.0 );
        g[i] = f(x[i]) * step_size * step_size;
    }

    double btemp = b[0];
    v[0] = g[0] / btemp;
    for ( int i = 1; i < n - 1; ++i ) {
        temp[i] = c[i-1] / btemp;
        btemp = b[i] - a[i] * temp[i];
        v[i] = ( g[i] - a[i] * v[i-1] ) / btemp;
    }

    for ( int i = n - 2; i >= 0; --i ) {
        v[i] -= temp[i+1] * v[i + 1];
    }

    for ( int i = 0; i < n; ++i ) {
        std::cout << v[i] << std::endl;
    }

    return 0;
}

double f(double x)
{
    return 100 * exp(-10*x);
}
