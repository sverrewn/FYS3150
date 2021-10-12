#include <armadillo>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "penning_trap.hpp"
#include "particle.hpp"


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

    return {x, y, z};
}


arma::vec PenningTrap::external_B_field(arma::vec r)
{
    return {0, 0, B};
}


arma::vec PenningTrap::particle_E_field(int i, int j)
{   
    arma::vec res = particles[i].position - particles[j].position;
    double qj = particles[j].charge;

    return ( ke * qj * res / std::pow(arma::norm(res), 3) );
}

    
arma::vec PenningTrap::total_particle_E_field(int i)
{
    arma::vec tot_field(3);

    for ( int j = 0; j < i; ++j ) {
        tot_field += particle_E_field(i, j);
    }

    for ( long unsigned int j = i + 1; j < particles.size(); ++j ) {
        tot_field += particle_E_field(i, j);
    }

    return tot_field;
}


arma::vec PenningTrap::total_force_E_fields(int i)
{   
    Particle p = particles[i];
    double q = p.charge;
    arma::vec r = p.position;

    arma::vec tot_force = q * total_particle_E_field(i) + q * external_E_field(r);

    return tot_force;
}


arma::vec PenningTrap::total_force(int i)
{   
    Particle p = particles[i];
    double q = p.charge;

    arma::vec r = p.position;
    arma::vec v = particles[i].velocity;

    arma::vec B = external_B_field(r);


    arma::vec F = total_force_E_fields(i) + arma::cross(q*v, B);

    return F;
}


void PenningTrap::evolve_RK4(double dt)
{   
    arma::vec k1v, k1r, k2v, k2r, k3v, k3r, k4v, k4r;
    arma::vec r_old, r_new, v_old, v_new;
    arma::vec F_i;

    double m;

    for ( int i = 0; i < particles.size(); ++i ) {

        Particle& p = particles[i];
        F_i = total_force(i);

        m = p.mass;
        r_old = p.position;
        v_old = p.velocity;
        
        k1v = dt * F_i / m;
        k1r = dt * v_old;

        r_new = r_old + k1r / 2;
        v_new = v_old + k1v / 2;

        p.position = r_new;
        p.velocity = v_new;


        F_i = total_force(i);

        k2v = dt * F_i / m;
        k2r = dt * v_new;
        
        r_new = r_old + k2r / 2;
        v_new = v_old + k2v / 2;

        p.position = r_new;
        p.velocity = v_new;

        F_i = total_force(i);

        k3v = dt * F_i / m;
        k3r = dt * v_new;

        r_new = r_old + k3r / 2;
        v_new = v_old + k3v / 2;   
        
        p.position = r_new;
        p.velocity = v_new;

        F_i = total_force(i);

        k4v = dt * F_i / m;
        k4r = dt * v_new;
        
        p.position = r_old + (k1r + 2 * k2r + 2 * k3r + k4r) * 1/6;
        p.velocity = v_old + (k1v + 2 * k2v + 2 * k3v + k4v) * 1/6;
    }


    return;
}


void PenningTrap::evolve_euler_cromer(double dt)
{
    arma::vec F_i;
    double m;

    for ( int i = 0; i < particles.size(); ++i ) {
        Particle& p = particles[i];
        F_i = total_force(i);
        m = p.mass;
        p.velocity = p.velocity + dt * F_i/m;
        p.position = p.position + dt * p.velocity;
    }

    return;
}

void PenningTrap::print_particles() {
    for ( int i = 0; i < particles.size(); ++i ) {
        std::cout << "Particle " << i << std::endl;
        std::cout << "r: " << particles[i].position << std::endl;
        std::cout << "v: " << particles[i].velocity <<std::endl;
    }
}

void PenningTrap::write_particles(std::ofstream& file) {
    for ( int i = 0; i < particles.size(); ++i ) {
        Particle& p = particles[i];
        file << "particle " << i << std::endl;
        file << "r: " << p.position(0) << "," << p.position(1) << "," << p.position(2) << std::endl;
        file << "v: " << p.velocity(0) << "," << p.velocity(1) << "," << p.velocity(2) << std::endl;
    }
    file << std::endl;

}