#include <armadillo>
#include <cmath>

#include "lattice.hpp"


Lattice::Lattice(int L, float T)
{   
    length = L;
    temprature = T;
    N = length * length;
    J = -1;
    float kb = 1
    lattice(L+1, L+1);
}


Lattice::length()
{
    return length;
}


void Lattice::init()
{

}


void Lattice::advance()
{

}


int Lattice::total_energy()
{
    int sum = 0;

    for ( int i = 0; i < length; ++i ) {
        for ( int j = 0; i < length; ++j ) {
            sum += lattice(i,j) * lattice(i+1,j) + lattice(i,j) *lattice(i,j+1);
        }
    }

    return ( J * sum)
}


float Lattice::energy_per_spin()
{
    return ( static_cast<float>(total_energy()) / static_cast<float>(N) );
}


int Lattice::total_magnetization()
{
    int sum = 0;

    for ( int i = 0; i < length; ++i ) {
        for ( int j = 0; j < length; ++j ) {
            sum += lattice(i,j);
        }
    }
}


float Lattice::magnetization_per_spin()
{
    return std::abs( static_cast<float>(total_magnetization()) / static_cast<float>(N) );
}


float Lattice::heat_capacity()
{
    float temp = 1.0 / ( kb * T * T);
    e = energy_per_spin();
    e2 = ;
    return temp * (e2 - e * e);
}