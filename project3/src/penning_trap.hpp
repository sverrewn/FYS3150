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
    std::vector<Particle> particles;

    arma::vec external_E_field(arma::vec r);
    arma::vec external_B_field(arma::vec r);
    arma::vec particle_E_field(int i, int j);
    arma::vec total_particle_E_field(int i, bool particle_interact);
    arma::vec total_force_E_fields(int i, bool particle_interact);

public:
    PenningTrap();
    PenningTrap(double B, double V, double d);

    void add_particle(Particle p);
    void evolve_RK4(double dt, bool particle_interact);
    void evolve_euler_cromer(double dt, bool particle_interact);

    arma::vec total_force(int i, bool particle_interact);
    
    void print_particles();
    void write_particles(std::ofstream& file);
    std::vector<Particle> get_all_particles();
};

#endif