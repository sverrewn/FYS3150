#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

#include "lattice.hpp"


Lattice::Lattice(int L, float T, bool oredered)
    : lattice(L, L) // ugly c++ initialise before the constructor runs syntax. Armadillo is difficult :(
{   
    length = L;
    N = length * length;
    ordered = oredered;
    temperature = T;
    E = M = 0.;
    init(ordered);
}


// Wrap around for periodic boundary conditions
inline int Lattice::periodic_idx(int i)
{
    return ( (i + length) % length );
}


// Fill the lattice with random spins
void Lattice::init(bool ordered)
{   
    if ( ordered) {
        for ( int i = 0; i < length; ++i ) {
            for ( int j = 0; j < length; ++j ) {
                lattice(i,j) = 1;
                M += (lattice(i,j));
            }
        }
    }
    else {
        for ( int i = 0; i < length; ++i ) {
            for ( int j = 0; j < length; ++j ) {
                int num = arma::randi(arma::distr_param(0,1));
                if ( num < 1) {
                    lattice(i,j) = -1; 
                }
                else {
                    lattice(i,j) = 1;
                }
                M += (lattice(i,j));
            }
        }
    }

    for ( int i = 0; i < length; ++i ) {
        for ( int j = 0; j < length; ++j ) {
            E -= lattice.at(i,j) * (
                    lattice.at(i, periodic_idx(j + 1)) +
                    lattice.at(periodic_idx(i + 1), j)
                );
        }
    }

    e_look = std::vector<double>(17, 0.0);
    for ( int i = -8; i <= 8; i+=4 ) {
        e_look[i + 8] = exp(-i/temperature);
    }

    average = std::vector<double>(5, 0.);
}


void Lattice::metropolis()
{   
    for ( int i = 0; i < N; ++i ) {
        int x = arma::randi(arma::distr_param(0, length - 1));
        int y = arma::randi(arma::distr_param(0, length - 1));
        double deltaE = 2 * lattice.at(x,y) * (
            lattice.at(x, periodic_idx(y + 1)) +
            lattice.at(x, periodic_idx(y - 1)) +
            lattice.at(periodic_idx(x + 1), y) +
            lattice.at(periodic_idx(x - 1), y)
            );

        if ( arma::randu() <= e_look[deltaE + 8]) {
            lattice.at(x,y) *= -1;
            E += deltaE;
            M += 2 * lattice.at(x,y);
        }
    }
}


// advance n Monte Carlo Cycles
void Lattice::MCcycle(unsigned int n)
{
    for ( unsigned int i = 0; i < n; ++i ) {
        metropolis();
        average[0] += E; average[1] += E * E;
        average[2] += M; average[3] += M * M;
        average[4] += std::fabs(M);
    }

    write_results(n);
}


void Lattice::write_results(int cycles)
{   
    std::cout << "Cycles: " << cycles << " | T: " << temperature << " | Ord: " << ordered <<std::endl;
    std::cout << "E current: " << E << std::endl;
    std::cout << "M current: " << M << std::endl;
    std::cout << "<E>: "   << average[0]/cycles << std::endl;
    std::cout << "<E^2>: " << average[1]/cycles << std::endl;
    std::cout << "<M>: "   << average[2]/cycles << std::endl;
    std::cout << "<M^2>: " << average[3]/cycles << std::endl;
    std::cout << "<|M|>: " << average[4]/cycles << std::endl;
    float temp = 1.0 / ( N * temperature * temperature);
    std::cout << "Cv: "    << temp * (average[1]/cycles - average[0]/cycles * average[0]/cycles ) << std::endl;
    temp = 1.0 / (N * temperature);
    std::cout << "X: "     << temp * (average[3]/cycles - average[2]/cycles * average[2]/cycles ) << std::endl;
    std::cout << std::endl;
}


/* ############################### */
/*          Garbage heap           */
/* ############################### */

/*
// Get total energy of the Lattice
int Lattice::total_energy()
{
    int sum = 0;

    for ( int i = 0; i < length; ++i ) {
        for ( int j = 0; i < length; ++j ) {
            sum -= lattice(i,j) * lattice(i+1,j) + lattice(i,j) *lattice(i,j+1);
        }
    }

    return sum
}


// Get average energy per spin in the Lattice
float Lattice::energy_per_spin()
{
    return ( static_cast<float>(total_energy()) / static_cast<float>(N) );
}


// Get total magnetization of the Lattice
int Lattice::total_magnetization()
{
    int sum = 0;

    for ( int i = 0; i < length; ++i ) {
        for ( int j = 0; j < length; ++j ) {
            sum += lattice(i,j);
        }
    }
}


// Get average magnetization per spin in the Lattice
float Lattice::magnetization_per_spin()
{
    return std::abs( static_cast<float>(total_magnetization()) / static_cast<float>(N) );
}

*/