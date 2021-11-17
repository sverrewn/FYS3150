#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <armadillo>
#include <vector>


class Lattice {
private:
    int length;
    int temperature;
    int N;
    double E, M; // energy, magnetization
    std::vector e_look, average;
    arma::Mat<short> lattice;
public:
    Lattice(int L, float T);
    void init();
    int periodic_idx(int i);
    void metropolis();
    void advance(int n);
};


#endif