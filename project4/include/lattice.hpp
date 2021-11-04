#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <armadillo>

class Lattice {
private:
    int length;
    arma::mat<double> lattice;
public:
    Lattice(int L);
    int length();
};

#endif