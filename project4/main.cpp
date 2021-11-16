#include<iostream>
#include<armadillo>


std::vector<int> Delta_E_values(){
    // List possible Delta E values
    // Do not touch
    return {-8, -4, 0, 4, 8};
}

std::vector<double> calculate_possible_boltzman_factors(std::vector<int> &possible_Delta_E_values, double &beta){
    // Find possible values for the Boltzman factor
    std::vector<double> possible_Boltzman_factors;
    for (int i=0; i<possible_Delta_E_values.size(); i++){
        possible_Boltzman_factors.push_back(exp(-beta*(-possible_Delta_E_values[i])));
    }
    return possible_Boltzman_factors;
}

int main(){
    double beta = 1;
    std::vector<int> possible_Delta_E_values = Delta_E_values();
    std::vector<double> possible_Boltzman_factors = calculate_possible_boltzman_factors(possible_Delta_E_values, beta);


    return 0;
}