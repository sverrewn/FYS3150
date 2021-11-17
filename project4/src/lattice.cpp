#include <armadillo>
#include <cmath>
#include <random>

#include "lattice.hpp"


Lattice::Lattice(int L, float T)
{   
    length = L;
    temprature = T;
    N = length * length;
    E = M = 0.;
    init();
}


// Wrap around for periodic boundary conditions
inline int Lattice::periodic_idx(int i)
{
    return ( (i + L) % L );
}


// Fill the lattice with random spins
void Lattice::init()
{   
    lattice = arma::mat(L, L);

    // set up random number generation
    unsigned int seed = 123456789;
    mt19937 generator;
    generator.seed(seed);

    // generate either 0 or 1. Will substitute 0 with -1
    uniform_int_distribution<int> init_val(0, 1)

    for ( int i = 0; i < L; ++i ) {
        for ( int j = 0; j < L; ++j ) {
            int num = init_val(generator);
            if ( num < 1) {
                lattice(i,j) = -1; 
            }
            else {
                lattice(i,j) = 1;
            }
            M += static_cast<double>(lattice(i,j));
        }
    }

    for ( int i = 0; i < L; ++i ) {
        for ( int j = 0; j < L; ++j ) {
            E -= static_cast<double>(
                lattice.at(i,j) * (
                    lattice.at(i, periodic_idx(j + 1)) +
                    lattice.at(periodic_idx(i + 1), j)
                    )
                );
        }
    }

    e_look = std::vector(17, 0);
    for ( int i = -8; i <= 8; i+=4 ) {
        e_look[i + 8] = exp(-i/temprature);
    }

    average = std::vector(5, 0);
}


void Lattice::metropolis()
{
    for ( int i = 0; i < N; ++i ) {
        int x = arma::randi(arma::distr_param(0, L - 1));
        int y = arma::randi(arma::distr_param(0, L - 1));
        int deltaE = 2 * lattice.at(x,y) * (
            lattice.at(i, periodic_idx(j + 1)) +
            lattice.at(i, periodic_idx(j - 1)) +
            lattice.at(periodic_idx(i + 1), j) +
            lattice.at(periodic_idx(i - 1), j)
            );

        if ( arma::randi(arma::distr_param(0,1)) <= e_look[deltaE + 8]) {
            lattice.at(x,y) *= -1;
            E += static_cast<double>(deltaE);
            M += 2 * lattice.at(x,y);
        }
    }
}


// advance n Monte Carlo Cycles
void Lattice::advance(int n)
{
    for ( int i = 0; i < n; ++i ) {
        metropolis();
        average[0] += E; average[1] += E * E;
        average[2] += M; average[3] += M * M;
        average[4] += std::fabs(M);
    }
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


// Get specific heat capacity normalized to number of spins
float Lattice::heat_capacity()
{
    float temp = 1.0 / ( kb * T * T);
    e = energy_per_spin();
    e2 = ;
    return temp * (e2 - e * e);
}
*/