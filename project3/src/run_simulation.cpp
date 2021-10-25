#include <armadillo>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <vector>
#include <tuple>

#include "penning_trap.hpp"
#include "particle.hpp"


void broad_freq_scan(PenningTrap&);
void narrow_freq_scan(PenningTrap&, double, double);


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

    broad_freq_scan(initial_pt);
    narrow_freq_scan(initial_pt, 0.34, 0.71);

    return 0;
}


void broad_freq_scan(PenningTrap& initial_pt)
{
    PenningTrap pt;
    std::vector<std::tuple<double, int>> particles_remaining;

    std::vector<double> amplitudes {0.1, 0.4, 0.7};

    double dt = 0.005;
    int us = 500; // time

    static int iteration = 0;
    int total_its = 2.3/0.02 * 3 + 1; // Approximate
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
        file  << std::setw(18) << std::setprecision(10) << std::scientific << omega
              << std::setw(18) << std::setprecision(10) << std::scientific << remaining 
              << std::endl;
    }

    return;
}


void narrow_freq_scan(PenningTrap& initial_pt, double start, double end)
{
    std::vector<std::tuple<double, double, int>> particles_remaining_no_interact;
    std::vector<std::tuple<double, double, int>> particles_remaining_interact;

    std::vector<double> amplitudes {0.1, 0.4, 0.7};

    double dt = 0.005;
    int freq_dt = 5;
    int us = 500; // time
    int init_i = start * 1000;
    int end_i = end * 1000;

    int iteration = 0;
    int total_its =  3 * (end - start) / 0.005 + 1; // Very approximate. Will miss by a few iterations
    for ( auto f : amplitudes) {

        #pragma omp parallel for
        for ( int i = init_i; i <= end_i; i += freq_dt ) {
            PenningTrap pt = initial_pt;
            double o = i / 1000.0;

            pt.coloumb_switch(false);
            for ( double t = 0; t <= us; t += dt ) {
                pt.fluctuate_E_field(f, o, t);
                pt.evolve_RK4(dt);
            }
            #pragma omp critical
            {
                particles_remaining_no_interact.push_back(std::make_tuple(f, o, pt.particles_trapped()));
                std::cout << "iteration " << ++iteration << "/" << total_its << std::endl;
            }
        }

    }
    std::ofstream file1;
    file1.open("data/narrow_freq_no_interact.dat");
    double freq; double omega; int remaining;
    for ( auto pair : particles_remaining_no_interact) {
        std::tie (freq, omega, remaining) = pair;
        file1 << std::setw(18) << std::setprecision(10) << std::scientific << freq
              << std::setw(18) << std::setprecision(10) << std::scientific << omega
              << std::setw(18) << std::setprecision(10) << std::scientific << remaining 
              << std::endl;
    }

    std::cout << "Starting round 2 with particle interaction" << std::endl;
    iteration = 0;
    for ( auto f : amplitudes) {

        #pragma omp parallel for
        for ( int i = init_i; i <= end_i; i += freq_dt ) {
            PenningTrap pt = initial_pt;
            double o = i / 1000.0;

            for ( double t = 0; t <= us; t += dt ) {
                pt.fluctuate_E_field(f, o, t);
                pt.evolve_RK4(dt);
            }
            #pragma omp critical
            {
                particles_remaining_interact.push_back(std::make_tuple(f, o, pt.particles_trapped()));
                std::cout << "iteration " << ++iteration << "/" << total_its << std::endl;
            }
        }

    }

    std::ofstream file2;
    file2.open("data/narrow_freq_interact.dat");

    for ( auto pair : particles_remaining_interact) {
        std::tie (freq, omega, remaining) = pair;
        file2 << std::setw(18) << std::setprecision(10) << std::scientific << freq
              << std::setw(18) << std::setprecision(10) << std::scientific << omega
              << std::setw(18) << std::setprecision(10) << std::scientific << remaining 
              << std::endl;
    }

    return;
}