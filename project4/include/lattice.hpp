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
    void init(bool ordered);
    int periodic_idx(int i);
    void metropolis();
    void write_results(int cycles, std::string base_name);
public:
    Lattice(int L, float T, bool ordered);
    void MCcycle(unsigned int n, std::string base_name);
    void MCcycle_no_write(unsigned int n);
    void MCcycle_n_samples_eps(unsigned int n, std::string fname);

};


#endif