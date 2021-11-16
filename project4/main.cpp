#include<iostream>
#include<armadillo>


std::vector<int> Delta_E_values(){
    // List possible Delta E values
    // Do not touch
    return {-8, -4, 0, 4, 8};
}

int main(){
    std::vector<int> possible_Delta_E_values = Delta_E_values();

    return 0;
}