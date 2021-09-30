#include "particle.hpp"

Particle::Particle(double c, double m, arma::vec p, arma::vec v) {
    charge = c;
    mass = m;
    position = p;
    velocity = v;
}