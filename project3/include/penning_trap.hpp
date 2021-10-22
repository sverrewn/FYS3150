#ifndef PENNING_TRAP_HPP
#define PENNING_TRAP_HPP

#include <armadillo>
#include <fstream>
#include <vector>
#include "particle.hpp"


class PenningTrap {
private:
    double B, V, V0;
    double d;
    double ke;
    std::vector<Particle> particles;

    bool coloumb_interact;

    // Private helpers
    arma::vec total_force_E_fields(int i);
    arma::vec total_particle_E_field(int i);
    arma::vec particle_E_field(int i, int j);
    arma::vec external_E_field(arma::vec r);
    arma::vec external_B_field(arma::vec r);

public:
    PenningTrap();
    PenningTrap(double V_in, double d_in);

    // Methods for modifying the trap
    void add_particle(Particle p);
    void fluctuate_E_field(double f, double o, double t);
    void coloumb_switch(bool value);
    void evolve_RK4(double dt);
    void evolve_euler_cromer(double dt);

    // Get the total force
    arma::vec total_force(int i);
    
    // Methods to extract information about the trap
    double get_d();
    void print_particles();
    void write_particles(std::ofstream& pos_file, std::ofstream& vel_file);
    int particles_trapped();
    std::vector<Particle> get_all_particles();
};

#endif