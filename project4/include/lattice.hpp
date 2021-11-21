#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <armadillo>
#include <string>
#include <vector>


class Lattice {
private:
    int length;
    int N;
    bool ordered;
    float temperature;
    double E, M; // energy, magnetization
    std::vector<double> e_look, average;
    arma::Mat<double> lattice;
public:
    Lattice(int L, float T, bool ordered);
    void init(bool ordered);
    int periodic_idx(int i);
    void metropolis();
    void MCcycle(unsigned int n, std::string base_name);
    void write_results(int cycles, std::string base_name);
};


#endif