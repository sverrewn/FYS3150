#include <armadillo>
#include <fstream>
#include <iostream>
#include <vector>

#include "penning_trap.hpp"
#include "particle.hpp"

int main()
{
    arma::arma_rng::set_seed(123456);
    // Verifying that RK4 and Euler-cromer are more or less in agreement
    
    PenningTrap pt_rk4 = PenningTrap();
    PenningTrap pt_ec = PenningTrap();

    pt_rk4.add_particle(Particle({2.7444e+01, -1.2353e+02, 7.4481e+01},
                                 {53.6670, -62.4112, 34.6811}));
    pt_rk4.add_particle(Particle({-58.4062, 23.6787, -46.4470},
                                 {-1.6068e+01, -1.2574e+02, 3.1104e+01}));
 
    pt_ec.add_particle(Particle({2.7444e+01, -1.2353e+02, 7.4481e+01},
                                 {53.6670, -62.4112, 34.6811}));
    pt_ec.add_particle(Particle({-58.4062, 23.6787, -46.4470},
                                {-1.6068e+01, -1.2574e+02, 3.1104e+01}));
    double dt = 0.01;
    std::ofstream file_rk4, file_ec;
    file_rk4.open("rk4_res.txt");
    file_ec.open("ec_res.txt");

    for ( int i = 0; i < 5; ++i ) {
        std::cout << "iteration " << i << std::endl;
        pt_rk4.evolve_RK4(dt);
        pt_ec.evolve_euler_cromer(dt);

        pt_rk4.write_particles(file_rk4);
        pt_ec.write_particles(file_ec);
    }

    return 0;    
}

