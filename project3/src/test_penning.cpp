#include <armadillo>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "penning_trap.hpp"
#include "particle.hpp"


void test_one_particle();
void test_two_particles();
void compare_RK4_EC();
void gen_analytical();


int main()
{   
    test_one_particle();
    test_two_particles();
    //compare_RK4_EC();
    gen_analytical();
    return 0;    
}


void test_one_particle()
{
    arma::arma_rng::set_seed(123456);
    PenningTrap pt = PenningTrap();
    
    pt.add_particle(Particle(arma::vec(3).randn() * 0.1 * 1000, arma::vec(3).randn() * 0.1 * 1000));

    std::ofstream file;
    file.open("data/single_particle_100us.txt");

    double dt = 0.001;
    int n_steps = 100 * 1000;
    pt.write_particles(file);
    for ( int i = 0; i < n_steps; ++i ) {
        pt.evolve_RK4(dt, true);
        pt.write_particles(file);
    }

    return;
}


void test_two_particles()
{
    arma::arma_rng::set_seed(3123132132);
    arma::vec r1 = arma::vec(3).randn() * 0.1 * 1000;
    arma::vec v1 = arma::vec(3).randn() * 0.1 * 1000;
    arma::vec r2 = arma::vec(3).randn() * 0.1 * 1000;
    arma::vec v2 = arma::vec(3).randn() * 0.1 * 1000;
    
    PenningTrap pt_interact = PenningTrap();
    PenningTrap pt_no_interact = PenningTrap();

    pt_interact.add_particle({r1,v1});
    pt_interact.add_particle({r2, v2});

    pt_no_interact.add_particle({r1, v1});
    pt_no_interact.add_particle({r2, v2});

    std::ofstream file1, file2;
    file1.open("data/two_particles_interaction.txt");
    file2.open("data/two_particles_no_interaction.txt");

    double dt = 0.001;
    int n_steps = 100 * 1000;    
    for ( int i = 0; i < n_steps; ++i ) {
        pt_interact.evolve_RK4(dt, true);
        pt_no_interact.evolve_RK4(dt, false);

        pt_interact.write_particles(file1);
        pt_no_interact.write_particles(file2);
    }

    return;
}


void compare_RK4_EC()
{   
    arma::arma_rng::set_seed(123456);
    arma::vec r1 = arma::vec(3).randn() * 0.1 * 1000;
    arma::vec v1 = arma::vec(3).randn() * 0.1 * 1000;

    std::vector<float> dt = {1, 0.1, 0.01, 0.001, 0.0005};

    // simulation time in us
    int time = 300;

    for ( auto h : dt ) {
        std::vector<Particle> p_rk4, p_ec;
        arma::vec avg_rel_err;

        int dt_per_time_unit = 1 / h;
        int n_steps = time * dt_per_time_unit;

        PenningTrap pt_rk4 = PenningTrap();
        PenningTrap pt_ec = PenningTrap();

        pt_rk4.add_particle({r1, v1}); 
        pt_ec.add_particle({r1, v1});

        // get nice numbers in the text file name
        std::string num = std::to_string(h);
        num.erase(num.find_last_not_of('0') + 1, std::string::npos);
        num.erase(num.find_last_not_of('.') + 1, std::string::npos);

        std::string fname1 = "data/one_particle_rk4_dt" + num + ".txt";
        std::string fname2 = "data/one_particle__ec_dt" + num + ".txt";
        std::ofstream file1, file2;
        file1.open(fname1);
        file2.open(fname2);

        pt_rk4.write_particles(file1);
        pt_ec.write_particles(file2);

        for ( int i = 0; i < n_steps; ++i ) {
            pt_rk4.evolve_RK4(h, false);
            pt_ec.evolve_euler_cromer(h, false);

            pt_rk4.write_particles(file1);
            pt_ec.write_particles(file2);

        }
    }

    return;
}


void gen_analytical()
{
    
}