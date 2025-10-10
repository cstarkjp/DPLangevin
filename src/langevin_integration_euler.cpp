/**
 * @file langevin_integration_euler.cpp
 * @brief Methods to carry out integration by explicit-Euler time-stepping.
 */ 

#include "general_core.hpp"

//! Perform explicit-Euler then stochastic integration steps, and update grid
void BaseLangevin::integrate_euler(rng_t &rng)
{
    euler_and_stochastic(density_grid_aux_old, rng);
    // Update density field grid with result of integration
    density_grid.swap(density_grid_aux_old); 
}

//! Perform an explicit-Euler integration step then a stochastic integration 
//! step for all cells across the grid.
void BaseLangevin::euler_and_stochastic(grid_t &density_grid_aux, rng_t &rng)
{
    mean_density = 0.0;
    for (auto i=0; i<n_cells; i++)
    {
        double f = nonlinear_rhs(i, density_grid);
        density_grid_aux[i] = density_grid[i] + f*dt;

        #if !APPROXIMATE_POISSON_DISTBN
        poisson_rng = poisson_dist_t(lambda_product * density_grid_aux[i]);
        gamma_rng = gamma_dist_t(poisson_rng(rng), 1.0);
        #else
        double mu = lambda_product * density_grid_aux[i];
        if (mu > MU_THRESHOLD)
        {
            normal = normal_dist_t(mu, sqrt(mu));
            gamma_rng = gamma_dist_t(normal(rng), 1.0);
        }
        else
        {
            poisson_rng = poisson_dist_t(lambda_product * density_grid_aux[i]);
            gamma_rng = gamma_dist_t(poisson_rng(rng), 1.0);
        }
        #endif
        density_grid_aux[i]= gamma_rng(rng)/lambda;
        mean_density += density_grid_aux[i];
    }    
    mean_density /= static_cast<double>(n_cells);   
}