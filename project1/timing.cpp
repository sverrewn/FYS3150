#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "algorithms.hpp"

double f(double);


int main()
{
    const std::vector<int> n_vals = {100, 1'000, 10'000, 100'000, 1'000'000};
    std::vector<std::tuple<int, double>> gen_results;
    std::vector<std::tuple<int, double>> spec_results;

    std::cout << "Testing general_tridiag\n" << std::endl;

    for ( auto n: n_vals ) {
        double result [3];

        std::cout << "Testing n=" << n << std::endl;

        for ( int i = 0; i < 3; ++i ) {
            std::vector<double>  a(n, -1), b(n, 2), c(n, -1), x, v, g, temp;
            x.reserve(n); v.reserve(n); g.reserve(n), temp.reserve(n);

            const double step_size = 1.0 / (n - 1);
            for ( int i = 0; i < n; ++i ) {
                x[i] = i / ( n - 1.0 );
                g[i] = f(x[i]) * step_size * step_size;
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            general_tridiag(n, a, b, c, v, g, temp);
            auto t2 = std::chrono::high_resolution_clock::now();

            double duration_seconds = std::chrono::duration<double>(t2-t1).count();

            result[i] = duration_seconds;
        }

        double average = ( result[0] + result[1] + result[2] ) / 3;
        gen_results.push_back(std::make_tuple(n, average));
    }
    std::cout << "\n-----------------------\n" << std::endl;
    std::cout << "Testing special_tridiag\n" << std::endl;

    for ( auto n: n_vals ) {
        double result [3];

        std::cout << "Testing n=" << n << std::endl;

        for ( int i = 0; i < 3; ++i ) {
            std::vector<double>  a(n, -1), b(n, 2), c(n, -1), x, v, g, temp;
            x.reserve(n); v.reserve(n); g.reserve(n), temp.reserve(n);

            const double step_size = 1.0 / (n - 1);
            for ( int i = 0; i < n; ++i ) {
                x[i] = i / ( n - 1.0 );
                g[i] = f(x[i]) * step_size * step_size;
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            general_tridiag(n, a, b, c, v, g, temp);
            auto t2 = std::chrono::high_resolution_clock::now();

            double duration_seconds = std::chrono::duration<double>(t2-t1).count();

            result[i] = duration_seconds;
        }

        double average = ( result[0] + result[1] + result[2] ) / 3;
        spec_results.push_back(std::make_tuple(n, average));
    }

    std::ofstream file;
    file.open("exec_time_cmp.dat");

    file << gen_results.size() << std::endl;

    for (auto n: gen_results) {
        file << std::get<0>(n) << "," << std::get<1>(n) << std::endl;
    }

    std::cout << "\n special times" << std::endl;
    for (auto n: spec_results) {
        file << std::get<0>(n) << "," << std::get<1>(n) << std::endl;
    }
}

double f(double x)
{
    return 100 * exp(-10*x);
}
