#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <armadillo>

class Particle {
private:
    double charge;
    double mass;
    arma::vec position;
    arma::vec velocity;
public:
    Particle(double c, double m, arma::vec p, arma::vec v);
    Particle(arma::vec p, arma::vec v);

    friend class PenningTrap;
};

#endif