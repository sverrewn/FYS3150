#include <cmath>
#include <vector>


double u(double);


void poisson_exact(const int n, const std::vector<double>& x_vals, std::vector<double>& result)
{
    for ( int i = 0; i < n; ++i ) {
        result.push_back(u(x_vals[i]));
    }
}


void general_tridiag(const int n, std::vector<double>& a, std::vector<double>& b,
                         std::vector<double>& c, std::vector<double>& v,
                         std::vector<double>& g, std::vector<double>& temp)
{
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
}

void special_tridiag(const int n, std::vector<double>& v,std::vector<double>& g,
                         std::vector<double>& temp)
{
    double btemp = 2;
    v[0] = g[0] / btemp;
    for ( int i = 1; i < n - 1; ++i ) {
        temp[i] = -1 / btemp;
        btemp = 2 + temp[i];
        v[i] = ( g[i] + v[i-1] ) / btemp;
    }

    for ( int i = n - 2; i >= 0; --i ) {
        v[i] -= temp[i+1] * v[i + 1];
    }
}

// Helper for poisson_exact
double u(double x)
{
    static const double u = exp(-10);
    return ( 1 - (1 - u) * x - exp(-10 * x) );
}
