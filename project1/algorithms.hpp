#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

void poisson_exact(const int n, const std::vector<double>& x_vals, std::vector<double>& result);

void general_tridiag(const int n, std::vector<double>& a, std::vector<double>& b,
                         std::vector<double>& c, std::vector<double>& v,
                         std::vector<double>& g, std::vector<double>& temp);

void special_tridiag(const int n, std::vector<double>& v,std::vector<double>& g,
                         std::vector<double>& temp);

#endif