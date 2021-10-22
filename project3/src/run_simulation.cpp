#include <armadillo>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <tuple>

#include "penning_trap.hpp"
#include "particle.hpp"


int main()
{
    PenningTrap initial_pt = PenningTrap(0.0025, 0.05);
    
    arma::arma_rng::set_seed(749257);
    int n = 100;

    for ( int i = 0; i < n; ++i ) {
        arma::vec r = arma::vec(3).randn() * 0.1 * initial_pt.get_d();
        arma::vec v = arma::vec(3).randn() * 0.1 * initial_pt.get_d();

        initial_pt.add_particle({r,  v});
    }


    PenningTrap pt;
    std::vector<std::tuple<double, int>> particles_remaining;

    std::vector<double> amplitudes {0.1, 0.4, 0.7};

    double dt = 0.005;
    int us = 500; // time

    static int iteration = 0;
    int total_its = 2.3/0.02 * 3 + 1;
    for ( auto f : amplitudes) {

        for ( double o = 0.2; o <= 2.5; o += 0.02 ) {
            pt = initial_pt;
            pt.coloumb_switch(false);

            for ( double t = 0; t <= us; t += dt ) {
                pt.fluctuate_E_field(f, o, t);
                pt.evolve_RK4(dt);
            }
            particles_remaining.push_back(std::make_tuple(o, pt.particles_trapped()));
            std::cout << "iteration " << ++iteration << "/" << total_its << std::endl;
        }

    }

    std::ofstream file;
    file.open("data/particles_remaining_for_f0.1_0.4_0.7.dat");

    double omega; int remaining;
    for ( auto pair : particles_remaining) {
        std::tie (omega, remaining) = pair;
        file  << std::setw(6) << std::setprecision(10) << std::scientific << omega
              << std::setw(6) << std::setprecision(10) << std::scientific << remaining 
              << std::endl;
    }
}