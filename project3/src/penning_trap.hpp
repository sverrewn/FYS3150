#ifndef PENNING_TRAP_HPP
#define PENNING_TRAP_HPP

#include <vector>
#include "particle.hpp"

class PenningTrap {
private:
    double B;
    double V;
    double d;
    std::vector<Particle> particles;
public:

    PenningTrap();
    PenningTrap(double B, double V, double d);

    void add_particle(Particle p);

    arma::vec external_E_field(arma::vec r);
    arma::vec external_B_field(arma::vec r);

    arma::vec force_particle(int i, int j);

    arma::vec total_force_external(int i);
    arma::vec total_force_particles(int i);

    arma::vec total_force(int i);

    void evolve_RK4(double dt);
    void evolve_forward_euler(double dt);

    virtual ~PenningTrap();
}

#endif