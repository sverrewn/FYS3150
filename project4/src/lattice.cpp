#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "lattice.hpp"


// Takes in size and temperature. ordered is whether or not the initial state is ordered (true), or random (false)
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


// Set up the Lattice
void Lattice::init(bool ordered)
{   
    // ordered spins
    if ( ordered) {
        for ( int i = 0; i < length; ++i ) {
            for ( int j = 0; j < length; ++j ) {
                lattice(i,j) = 1;
                M += (lattice(i,j));
            }
        }
    }
    else { // unordered spins
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

    // Calculate initial E
    for ( int i = 0; i < length; ++i ) {
        for ( int j = 0; j < length; ++j ) {
            E -= lattice.at(i,j) * (
                    lattice.at(i, periodic_idx(j + 1)) +
                    lattice.at(periodic_idx(i + 1), j)
                );
        }
    }

    // set up loopup table for probability of a Delta E
    e_look = std::vector<double>(17, 0.0);
    for ( int i = -8; i <= 8; i+=4 ) {
        e_look[i + 8] = exp(-i/temperature);
    }

    // vector to store results
    average = std::vector<double>(5, 0.);

    return;
}

// Advances the system on cycle, calculating N new possible states, and accepting them according to the metropolis rule
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


// Advance the system n MC cycles and write the results to a file
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


// Samples the eps value for n MC cycles
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

// Writes the results to a file
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
