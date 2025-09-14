// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

// Make the density visible
dbl_vec_t LangevinBase::get_density() 
{
    return cell_density;
}

// Make the density visible
double LangevinBase::get_mean_density() 
{
    return mean_density;
}

// Check the mean of the Poisson distribution
double LangevinBase::get_poisson_mean()
{
    return lambda_product * mean_density;
}