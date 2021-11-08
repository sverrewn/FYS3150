#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <armadillo>

class Lattice {
private:
    int length;
    int temperature;
    int N;
    int J;
    float kb;
    arma::Mat<short> lattice;
public:
    Lattice(int L, float T,);
    int length();

    float energy_per_spin();
};

#endif