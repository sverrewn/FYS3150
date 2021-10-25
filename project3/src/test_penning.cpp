#include <armadillo>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "penning_trap.hpp"
#include "particle.hpp"


void test_one_particle_random_initial_values();
void test_one_particle(arma::vec init_pos, arma::vec init_vel, std::string file_pos, std::string file_vel);
void test_two_particles();
void compare_RK4_EC();
void gen_analytical();
void test_one_particle_given_initial_values();


int main()
{   
    test_one_particle_random_initial_values();
    test_one_particle_given_initial_values();
    test_two_particles();
    compare_RK4_EC();
    
    return 0;    
}


void test_one_particle_random_initial_values()
{
    arma::arma_rng::set_seed(123456);
    arma::vec init_pos = arma::vec(3).randn() * 0.1 * 200, init_vel = arma::vec(3).randn() * 0.1 * 200;
    std::string pos_file = "data/single_particle_100us_pos.txt", vel_file = "data/single_particle_100us_vel.txt";
    test_one_particle(init_pos, init_vel, pos_file, vel_file);

    return;
}

void test_one_particle(arma::vec init_pos, arma::vec init_vel, std::string file_pos, std::string file_vel)
{
    PenningTrap pt = PenningTrap();
    
    pt.add_particle(Particle(init_pos, init_vel));

    std::ofstream pos_file, vel_file;
    pos_file.open(file_pos);
    vel_file.open(file_vel);

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
    double x_0 = 500, z_0 = 300, v_0 = 110;
    arma::vec r1 = {x_0, 0, z_0}, v1 = {0, v_0, 0};

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


void test_one_particle_given_initial_values()
{
    double x_0 = 500, z_0 = 300, v_0 = 110;
    arma::vec init_pos = {x_0, 0, z_0}, init_vel = {0, v_0, 0};
    std::string pos_file = "data/analytical_pos_test.txt", vel_file = "data/analytical_vel_test.txt";
    test_one_particle(init_pos, init_vel, pos_file, vel_file);
}