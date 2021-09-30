#include "penning_trap.hpp"

PenningTrap::PenningTrap() {
    B = 1.0;
    V = 10.0;
    d = 10000;
}

PenningTrap::PenningTrap(double B, double V, double d) {
    B = B;
    V = V;
    d = d;
}