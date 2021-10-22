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
    compare_RK4_EC();
    //gen_analytical();
    return 0;    
}


void test_one_particle()
{
    arma::arma_rng::set_seed(123456);
    PenningTrap pt = PenningTrap();
    
    pt.add_particle(Particle(arma::vec(3).randn() * 0.1 * 200, arma::vec(3).randn() * 0.1 * 200));

    std::ofstream pos_file, vel_file;
    pos_file.open("data/single_particle_100us_pos.txt");
    vel_file.open("data/single_particle_100us_vel.txt");

    double dt = 1e-3;
    double time = 100; 
    int n_steps = time * 1/dt;

    pt.write_particles(pos_file, vel_file);
    for ( int i = 0; i < n_steps; ++i ) {
        pt.evolve_RK4(dt);
        pt.write_particles(pos_file, vel_file);
    }

    return;
}


void test_two_particles()
{
    arma::arma_rng::set_seed(3123132132);
    arma::vec r1 = arma::vec(3).randn() * 50;
    arma::vec v1 = arma::vec(3).randn() * 50;
    arma::vec r2 = arma::vec(3).randn() * 50;
    arma::vec v2 = arma::vec(3).randn() * 50;
    
    PenningTrap pt_interact = PenningTrap();
    PenningTrap pt_no_interact = PenningTrap();

    pt_interact.add_particle({r1,v1});
    pt_interact.add_particle({r2, v2});

    pt_no_interact.add_particle({r1, v1});
    pt_no_interact.add_particle({r2, v2});
    pt_no_interact.coloumb_switch(false);

    std::ofstream pos_file1, vel_file1, pos_file2, vel_file2;
    pos_file1.open("data/two_particles_interaction_pos.txt");
    vel_file1.open("data/two_particles_interaction_vel.txt");
    pos_file2.open("data/two_particles_no_interaction_pos.txt");
    vel_file2.open("data/two_particles_no_interaction_vel.txt");

    double dt = 0.001;
    int n_steps = 100 * 1000;
    
    pt_interact.write_particles(pos_file1, vel_file1);
    pt_no_interact.write_particles(pos_file2, vel_file2);
    for ( int i = 0; i < n_steps; ++i ) {
        pt_interact.evolve_RK4(dt);
        pt_no_interact.evolve_RK4(dt);

        pt_interact.write_particles(pos_file1, vel_file1);
        pt_no_interact.write_particles(pos_file2, vel_file2);
    }

    return;
}


void compare_RK4_EC()
{   
    arma::arma_rng::set_seed(123456);
    arma::vec r1 = arma::vec(3).randn() * 0.1 * 1000;
    arma::vec v1 = arma::vec(3).randn() * 0.1 * 1000;

    std::vector<float> dt = {1, 1e-1, 1e-2, 1e-3, 5e-4};

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

        pt_rk4.coloumb_switch(false);
        pt_ec.coloumb_switch(false);

        // get nice numbers in the text file name
        std::string num = std::to_string(h);
        num.erase(num.find_last_not_of('0') + 1, std::string::npos);
        num.erase(num.find_last_not_of('.') + 1, std::string::npos);

        std::string fname1 = "data/one_particle_rk4_dt" + num + "_pos.txt";
        std::string fname2 = "data/one_particle_rk4_dt" + num + "_vel.txt";
        std::string fname3 = "data/one_particle__ec_dt" + num + "_pos.txt";
        std::string fname4 = "data/one_particle__ec_dt" + num + "_vel.txt";

        std::ofstream file1, file2, file3, file4;
        file1.open(fname1);
        file2.open(fname2);
        file3.open(fname3);
        file4.open(fname4);

        pt_rk4.write_particles(file1, file2);
        pt_ec.write_particles(file3, file4);

        for ( int i = 0; i < n_steps; ++i ) {
            pt_rk4.evolve_RK4(h);
            pt_ec.evolve_euler_cromer(h);

            pt_rk4.write_particles(file1, file2);
            pt_ec.write_particles(file3, file4);

        }
    }

    return;
}


void gen_analytical()
{
    return;    
}
