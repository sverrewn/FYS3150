#include <armadillo>
#include <fstream>
#include <iostream>
#include <vector>

#include "penning_trap.hpp"
#include "particle.hpp"


int main()
{
    PenningTrap initial_pt = PenningTrap(0.0025, 0.05);
    
    arma::arma_rng::set_seed(749257);
    int n = 100;

    for ( int i = 0; i < n; ++i ) {
        vec r = vec(3).randn() * 0.1 * initial_pt.get_d();
        vec v = vec(3).randn() * 0.1 * initial_pt.get_d();

        initial_pt.add_particle({r,  v})
    }


    PenningTrap pt;
    std::vec<double> amplitudes {0.1, 0.4, 0.7}

    std::string fname1, fname2;
    std::ofstream pos_file, vel_file;
    
    int us = 500; // time

    for ( auto f : amplitudes) {
        std::string num = std::to_string(f);
        num.erase(num.find_last_not_of('0') + 1, std::string::npos);
        num.erase(num.find_last_not_of('.') + 1, std::string::npos);
        
        fname1 = "f" + num + "_pos.txt";
        fname2 = "f" + num + "_vel.txt";
        pos_file1.open(fname1);
        vel_file1.open(fname2);
        
        pt = initial_pt;
        pt.coloumb_switch(false);


    }

}

