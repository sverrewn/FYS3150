#include <iostream>
#include <vector>

#include "penning_trap.hpp"
#include "particle.hpp"

int main()
{
    PenningTrap pt = PenningTrap();
    pt.add_particle({{0.1, 0.1, 0.05}, {1,1,1}});
    pt.add_particle({{0.05, 0.05, 0.05}, {1.1,1.1,1.1}});
    pt.evolve_RK4(0.001);
    return 0;    
}