// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "declarations.hpp"

void Dornic::construct_custom_network(const std::vector<int_vector> &network)
{
    neighbors = network;
}
