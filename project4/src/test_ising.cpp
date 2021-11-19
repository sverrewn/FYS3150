#include <iostream>
#include <vector>
#include <omp.h>
#include "lattice.hpp"


void verify();
void burn_in();
void find_approx_distr();

int main()
{
    //verify();

    //burn_in();

    find_approx_distr();


    return 0;
}


void verify()
{
    Lattice ising_2d = Lattice(2, 1, true);

    ising_2d.MCcycle(100 * 1000);
    
    return;
}


void burn_in()
{   
    int cycles[7] = {10, 100, 1000, 10000, 100000, 500000, 1000000};
    for ( int n : cycles) {
        Lattice ising_2d = Lattice(20, 1, true);
        ising_2d.MCcycle(n);

        ising_2d = Lattice(20, 1, false);
        ising_2d.MCcycle(n);

        ising_2d = Lattice(20, 2.4, true);
        ising_2d.MCcycle(n);

        ising_2d = Lattice(20, 2.4, false);
        ising_2d.MCcycle(n);
    }
}


void find_approx_distr()
{
    std::vector<double> temps;
    int n = 1000000;
    for ( int i = 100; i < 240; ++i ) {
        temps.push_back(i / 100.);
    }

    #pragma omp parallel for
    for ( int i = 0; i < temps.size(); ++i ) {
        Lattice ising_2d = Lattice(20, temps[i], false);
        ising_2d.MCcycle(n);
    }
    return;
}