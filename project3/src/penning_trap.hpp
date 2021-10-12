#ifndef PENNING_TRAP_HPP
#define PENNING_TRAP_HPP

#include <armadillo>
#include <fstream>
#include <vector>
#include "particle.hpp"

class PenningTrap {
private:
    double B;
    double V;
    double d;
    double vd2;
    double ke;

public:
    std::vector<Particle> particles;
    PenningTrap();
    
    PenningTrap(double B, double V, double d);

    void add_particle(Particle p);

    arma::vec external_E_field(arma::vec r);

    arma::vec external_B_field(arma::vec r);

    arma::vec particle_E_field(int i, int j);

    arma::vec total_particle_E_field(int i);
    
    arma::vec total_force_E_fields(int i);

    arma::vec total_force(int i);

    void evolve_RK4(double dt);
    
    void evolve_euler_cromer(double dt);

    void print_particles();

    void write_particles(std::ofstream& file);
};

#endif