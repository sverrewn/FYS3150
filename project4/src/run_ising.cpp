#include <iostream>
#include <omp.h>
#include <string>
#include <vector>


#include "lattice.hpp"


void verify(std::string);
void burn_in(std::string);
void test_speedup();
void find_approx_distr(int, std::string, std::string);
void find_phase(std::string);


int main()
{
    //verify("data/test/MC");

    //burn_in("data/burn_in/MC");

    //test_speedup();

    find_approx_distr(1000000, "data/approx_distr/MC1.txt", "data/approx_distr/MC2.4.txt");

    //find_phase("data/phase/");

    return 0;
}


void verify(std::string base_name)
{   
    std::cout << "Testing 2x2 Lattice" << "\n";
    int size = 6;
    int cycles[size] = {100, 1000, 10000, 100000, 1000000, 10000000};
    for ( int i = 0; i < size; ++i ) {
        std::cout << "\t" << cycles[i] << "...";
        Lattice ising_2d = Lattice(2, 1, false);
        ising_2d.MCcycle(cycles[i], base_name);
        std::cout << "done\n";
    }
    std::cout << std::endl;

    return;
}


void burn_in(std::string base_name)
{   
    std::cout << "Testing for burn_in\n";
    int cycles[5] = {100, 1000, 10000, 100000, 1000000};
    for ( int n : cycles) {
        std::cout << "\t" << n << "...";
        Lattice ising_2d = Lattice(20, 1, true);
        ising_2d.MCcycle(n, base_name);

        ising_2d = Lattice(20, 1, false);
        ising_2d.MCcycle(n, base_name);

        ising_2d = Lattice(20, 2.4, true);
        ising_2d.MCcycle(n, base_name);

        ising_2d = Lattice(20, 2.4, false);
        ising_2d.MCcycle(n, base_name);

        std::cout << "done\n";
    }
    std::cout << std::endl;
}


// MY MACHINE HAS 16 THREADS
void test_speedup()
{   
    std::cout << "Testing speed-up" << std::endl;
    double start, end;
    double results[6];

    
    for ( int i = 0; i < 3; ++ i ) {
        std::cout << "\ttesting ST " << i+1 << "/3" << std::endl;

        start = omp_get_wtime();
        for ( int i = 0; i < 16; ++i ) { 
            Lattice ising_2d = Lattice(20, 1, false);
            ising_2d.MCcycle_no_write(1000000);
        }
        end = omp_get_wtime();

        results[i] = end - start;
    }

    for ( int i = 3; i < 6; ++i ) {
        std::cout << "\tTesting MT " << i-2 << "/3" << std::endl;
        start = omp_get_wtime();

        #pragma omp parallel for
        for ( int i = 0; i < 16; ++i ) {
            Lattice ising_2d = Lattice(20, 1, false);
            ising_2d.MCcycle_no_write(1000000);
        }

        end = omp_get_wtime();

        results[i] = end - start;
    }

    double avg_st = ( results[0] + results[1] + results[2] ) / 3;
    double avg_mt = ( results[3] + results[4] + results[5] ) / 3;

    std::cout << "\nST speed: " << avg_st << "\n";
    std::cout << "MT speed: " << avg_mt << "\n";
    std::cout << "Speedup: " << (avg_mt/avg_st) * 100 << "\n" << std::endl;

}


void find_approx_distr(int n, std::string base_name1, std::string base_name2)
{   
    std::cout << "Sampling MCMC for T = 1...";

    Lattice ising_2d = Lattice(20, 1., false);
    ising_2d.MCcycle_n_samples_eps(n, base_name1);

    std::cout << "done" << std::endl;

    std::cout << "Sampling MCMC for T = 2.4";

    ising_2d = Lattice(20, 2.4, false);
    ising_2d.MCcycle_n_samples_eps(n, base_name2);

    std::cout << "done" << std::endl;

    return;
}


void find_phase(std::string base_name)
{
    std::vector<double> temps;
    int n = 200000;

    for ( int i = 2100; i < 2401; ++i ) {
        temps.push_back(i/1000.);
    }

    std::cout << "Starting L = 40" << std::endl;
    int approx_runs = 0;
    std::string name = base_name;
    name.append("L40");

    #pragma omp parallel
    {
    int tid = omp_get_thread_num();

    #pragma omp for
    for ( int i = 0; i < temps.size(); ++i ) {
        Lattice ising_2d = Lattice(40, temps[i], false);
        ising_2d.MCcycle(n, name);

        ++approx_runs;

        if ( tid == 0 ) {
            std::cout << "\t" << approx_runs << "/300 (approx)" << std::endl;
        }
    }
    }

    std::cout << "Starting L = 60" << std::endl;
    approx_runs = 0;
    name = base_name;
    name.append("L60");

    #pragma omp parallel
    { 
    int tid = omp_get_thread_num();

    #pragma omp for
    for ( int i = 0; i < temps.size(); ++i ) {
        Lattice ising_2d = Lattice(60, temps[i], false);
        ising_2d.MCcycle(n, name);

        ++approx_runs;

        if ( tid == 0 ) {
            std::cout << "\t" << approx_runs << "/300 (approx)" << std::endl;
        }

    }
    }

    std::cout << "Starting L = 80" << std::endl;
    approx_runs = 0;
    name = base_name;
    name.append("L80");

    #pragma omp parallel 
    {
    int tid = omp_get_thread_num();
    
    #pragma omp for
    for ( int i = 0; i < temps.size(); ++i ) {
        Lattice ising_2d = Lattice(80, temps[i], false);
        ising_2d.MCcycle(n, name);

        ++approx_runs;

        if ( tid == 0 ) {
            std::cout << "\t" << approx_runs << "/300 (approx)" << std::endl;
        }
    }
    }

    std::cout << "Starting L = 100" << std::endl;
    approx_runs = 0;
    name = base_name;
    name.append("L100");

    #pragma omp parallel 
    {
    int tid = omp_get_thread_num();

    #pragma omp for
    for ( int i = 0; i < temps.size(); ++i ) {
        Lattice ising_2d = Lattice(100, temps[i], false);
        ising_2d.MCcycle(n, name);

        ++approx_runs;

        if ( tid == 0 ) {
            std::cout << "\t" << approx_runs << "/300 (approx)" << std::endl;
        }
    }
    }
}