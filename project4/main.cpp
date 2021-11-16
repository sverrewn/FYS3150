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

void Metropolis(int L, long &idum, arma::mat &spins, std::vector<double> &possible_Boltzman_factors, double &E, double &M){
    // Loop over all spins
    for (int x = 0; x < L; x++){
        for (int y = 0; y < L; y++){
            std::cout << "Yay!";
            //Find random index values
            int i = random_number(idum)*L;
            int j = random_number(idum)*L;

            std::cout << "\n" << i << " " << j;
    
            int Delta_E = 2*spins.col(i)(j) * (spins.col(i)(j+1) + spins.col(i-1)(j) + spins.col(i)(j-1) + spins.col(i+1)(j));
            // Map Delta E to the corresponding index in possible_Boltzman_factors
            int index = Delta_E/4 + 2;
            double Boltzman_factor = possible_Boltzman_factors[index];

            if (random_number(idum) <= Boltzman_factor){
                spins.col(i)(j) *= -1;

                E += Delta_E;
                M += 2*spins.col(i)(j);
            }
        }
    }
}

int main(){
    double beta = 1;
    std::vector<int> possible_Delta_E_values = Delta_E_values();
    std::vector<double> possible_Boltzman_factors = calculate_possible_boltzman_factors(possible_Delta_E_values, beta);


    return 0;
}