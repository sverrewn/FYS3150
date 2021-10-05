#include <armadillo>
#include <cassert>
#include <vector>

#include "penning_trap.hpp"

PenningTrap::PenningTrap()
{
    B = 9.65 * 10; 
    V = std::pow(9.65, 8);
    d = 1000;

    assert ( d != 0);    
    vd2 = V / std::pow(d,2); // The ratio V/d^2

    ke = 1.38935333 * std::pow(10,5); // Coloumb's constant
}

PenningTrap::PenningTrap(double B0, double V0, double d)
{
    B = B0;
    V = V0;
    d = d;

    assert ( d != 0);    
    vd2 = V / std::pow(d,2);

    ke = 1.38935333 * std::pow(10,5); // Coloumb's constant
}

void PenningTrap::add_particle(Particle p)
{
    particles.push_back(p);
}

arma::vec PenningTrap::external_E_field(arma::vec r)
{
    double x = vd2 * r(0);
    double y = vd2 * r(1);
    double z = - 2 * vd2 * r(2);

    return {x,y,z};
}

arma::vec PenningTrap::external_B_field(arma::vec r)
{
    return {0,0, B};
}

arma::vec PenningTrap::force_particle(int i, int j)
{   
    arma::vec res = particles[i].position - particles[j].position;
    double qj = particles[j].charge;

    return ( ke * qj * res / std::pow(arma::norm(res), 3) );
}
    
arma::vec PenningTrap::total_force_particles(int i)
{
    arma::vec tot_field(3);

    for ( int j = 0; j < i; ++j ) {
        tot_field += force_particle(i, j);
    }

    for ( long unsigned int j = i + 1; j < particles.size(); ++j ) {
        tot_field += force_particle(i, j);
    }

    return tot_field;
}

arma::vec PenningTrap::total_force_external(int i)
{   
    Particle p = particles[i];
    double q = p.charge;
    arma::vec r = p.position;

    arma::vec tot_force = q * total_force_particles(i) + q * external_E_field(r);

    return tot_force;
}

arma::vec PenningTrap::total_force(int i)
{   
    Particle p = particles[i];
    double q = p.charge;

    arma::vec r = p.position;
    arma::vec v = particles[i].velocity;

    arma::vec B = external_B_field(r);
    arma::vec qvB = arma::vec( { q*v[1]*B[2], -B[2]*q*v[0], 0 } );
    arma::vec F = total_force_external(i) + qvB;
    
    return F;
}

void PenningTrap::evolve_RK4(double dt)
{
    return;
}
    
void PenningTrap::evolve_forward_euler(double dt)
{
    return;
}
