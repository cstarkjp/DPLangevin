// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

// Make the density visible
double LangevinBase::density() 
{
    return mean_density;
}

// Check the mean of the Poisson distribution
double LangevinBase::avg_poisson_mean()
{
    return lambda_product * mean_density;
}