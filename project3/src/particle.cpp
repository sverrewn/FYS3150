#include <cmath>
#include "particle.hpp"

Particle::Particle(double c, double m, arma::vec p, arma::vec v)
{
    charge = c;
    mass = m;
    position = p;
    velocity = v;
}

// Default Calcium ion
Particle::Particle(arma::vec p, arma::vec v)
{
    charge = 1;
    mass = 40.078;
    position = p;
    velocity = v;
}