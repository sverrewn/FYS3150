#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
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

    return;
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

    return;
}


// advance n Monte Carlo Cycles
void Lattice::MCcycle(unsigned int n, std::string base_name)
{
    for ( unsigned int i = 0; i < n; ++i ) {
        metropolis();
        average[0] += E; average[1] += E * E;
        average[2] += M; average[3] += M * M;
        average[4] += std::fabs(M);
    }

    write_results(n, base_name);

    return;
}


// This function is identical to MCcycle, but doesn't write to file.
// It exists solely to test speed without I/O interfering
void Lattice::MCcycle_no_write(unsigned int n)
{
    for ( unsigned int i = 0; i < n; ++i ) {
        metropolis();
        average[0] += E; average[1] += E * E;
        average[2] += M; average[3] += M * M;
        average[4] += std::fabs(M);
    }

    return;
}


void Lattice::MCcycle_n_samples_eps(unsigned int n, std::string fname)
{   
    std::ofstream file;
    file.open(fname);

    for ( unsigned int i = 1; i <= n; ++i ) {
        metropolis();

        file << std::setw(18) << std::setprecision(14) << E/N << std::endl;
    }

    return;
}


/*    Expected file layout
    -----------------------
    Ordered     (bool, 1/0)
    Cycles      (int)
    Temperature (double)
    Current E   (double)
    Current M   (double)
    <eps>       (double)
    <E^2>       (double)
    <M>         (double)
    <M^2>       (double)
    <|m|>       (double)
    Cv          (double)
    X           (double)
    -----------------------
*/

void Lattice::write_results(int cycles, std::string base_name)
{   
    std::string fname = base_name;
    fname.append("_T");
    fname.append(std::to_string(static_cast<int>(temperature*1000)));
    fname.append("_C");
    fname.append(std::to_string(cycles));
    fname.append(".txt");

    std::ofstream file;
    file.open(fname);
    
    double width = 18, prec = 14;    

    file << ordered << std::endl;
    file << std::setw(width) << std::setprecision(prec) << cycles << std::endl;
    file << std::setw(width) << std::setprecision(prec) << temperature << std::endl;
    file << std::setw(width) << std::setprecision(prec) << E << std::endl;
    file << std::setw(width) << std::setprecision(prec) << M << std::endl;
    file << std::setw(width) << std::setprecision(prec) << average[0]/(cycles * N) << std::endl;
    file << std::setw(width) << std::setprecision(prec) << average[1]/cycles << std::endl;
    file << std::setw(width) << std::setprecision(prec) << average[2]/cycles << std::endl;
    file << std::setw(width) << std::setprecision(prec) << average[3]/cycles << std::endl;
    file << std::setw(width) << std::setprecision(prec) << average[4]/(cycles * N) << std::endl;
    float temp = 1.0 / ( N * temperature * temperature);
    file << std::setw(width) << std::setprecision(prec) << temp * (average[1]/cycles - average[0]/cycles * average[0]/cycles ) << std::endl;
    temp = 1.0 / (N * temperature);
    file << std::setw(width) << std::setprecision(prec) << temp * (average[3]/cycles - average[4]/cycles * average[4]/cycles ) << std::endl;

    return;
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