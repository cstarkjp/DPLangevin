/**
 * @file langevin_integrate_euler.cpp
 * @brief Methods to carry out integration by explicit-Euler time-stepping.
 */ 

#include "langevin_types.hpp"
#include "langevin_base.hpp"

//! Perform explicit-Euler then stochastic integration steps, then update grid
void BaseLangevin::integrate_euler(rng_t& rng)
{
    mean_density = 0.0;
    for (auto i=0; i<n_cells; i++)
    {
        double f = nonlinear_rhs(i, density_grid);
        aux_grid1[i] = density_grid[i] + f*dt;
        poisson_rng = poisson_dist_t(lambda_product * aux_grid1[i]);
        gamma_rng = gamma_dist_t(poisson_rng(rng), 1.0);
        aux_grid1[i]= gamma_rng(rng)/lambda;
        mean_density += aux_grid1[i];
    }    
    mean_density /= static_cast<double>(n_cells);   
    // Update density field grid with result of integration
    density_grid.swap(aux_grid1); 
}