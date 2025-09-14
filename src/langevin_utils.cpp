// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "general_core.hpp"

// Make the density visible
dbl_vec_t Langevin::get_density() 
{
    return cell_density;
}

// Make the density visible
double Langevin::get_mean_density() 
{
    return mean_density;
}

// Check the mean of the Poisson distribution
double Langevin::get_poisson_mean()
{
    return lambda_product * mean_density;
}