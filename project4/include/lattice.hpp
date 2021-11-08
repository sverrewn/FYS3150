#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <armadillo>

class Lattice {
private:
    int length;
    arma::Mat<short> lattice;
public:
    Lattice(int L);
    int length();
};

#endif