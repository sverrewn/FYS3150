#include <iostream>
#include <vector>

#include "penning_trap.hpp"
#include "particle.hpp"

int main()
{
    PenningTrap pt = PenningTrap();
    pt.add_particle({{0.1, 0.1, 0.05}, {1,1,1}});
    pt.total_force(0);
    return 0;    
}