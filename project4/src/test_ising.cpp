#include <iostream>
#include <omp.h>
#include <string>
#include <vector>


#include "lattice.hpp"


void verify(std::string);
void burn_in(std::string);
void find_approx_distr(std::string);

int main()
{
    verify("data/test/MC");

    //burn_in();

   // find_approx_distr();


    return 0;
}


void verify(std::string base_name)
{
    Lattice ising_2d = Lattice(2, 1, true);

    ising_2d.MCcycle(100 * 1000, base_name);
    
    return;
}


void burn_in(std::string base_name)
{   
    int cycles[7] = {10, 100, 1000, 10000, 100000, 500000, 1000000};
    for ( int n : cycles) {
        Lattice ising_2d = Lattice(20, 1, true);
        ising_2d.MCcycle(n, base_name);

        ising_2d = Lattice(20, 1, false);
        ising_2d.MCcycle(n, base_name);

        ising_2d = Lattice(20, 2.4, true);
        ising_2d.MCcycle(n, base_name);

        ising_2d = Lattice(20, 2.4, false);
        ising_2d.MCcycle(n, base_name);
    }
}


void find_approx_distr(std::string base_name)
{
    std::vector<double> temps;
    int n = 1000000;
    for ( int i = 100; i < 240; ++i ) {
        temps.push_back(i / 100.);
    }

    #pragma omp parallel for
    for ( int i = 0; i < temps.size(); ++i ) {
        Lattice ising_2d = Lattice(20, temps[i], false);
        ising_2d.MCcycle(n, base_name);
    }
    return;
}