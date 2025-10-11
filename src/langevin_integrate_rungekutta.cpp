/**
 * @file langevin_integrate_rungekutta.cpp
 * @brief Methods to carry out 4th-order Runge-Kutta integration.
 */

#include "general_types.hpp"

//! Runge-Kutta integration of the nonlinear and diffusion terms 
//! in the Langevin equation.
//! Update of cells is done in the same loop as last Runge-Kutta step 
//! for efficiency.
void BaseLangevin::integrate_rungekutta(rng_t &rng)
{
    rk_f1(density_grid_aux_old, k1_grid);
    rk_f2f3(density_grid_aux_old, density_grid_aux_new, k2_grid, dt*0.5);
    // Swap contents is O(1), better than old = fast
    density_grid_aux_old.swap(density_grid_aux_new);            
    rk_f2f3(density_grid_aux_old, density_grid_aux_new, k3_grid, dt);
    density_grid_aux_old.swap(density_grid_aux_new);
    rk_f4_and_stochastic(density_grid_aux_old, k1_grid, k2_grid, k3_grid, rng);
}

//! Runge-Kutta first function
void BaseLangevin::rk_f1(grid_t &density_grid_aux, grid_t &k1_grid)
{
    for (auto i=0; i<n_cells; i++)
    {
        k1_grid[i] = nonlinear_rhs(i, density_grid);
        density_grid_aux[i] = density_grid[i] + 0.5*dt*k1_grid[i];
    }
}

//! Runge-Kutta second and third functions
void BaseLangevin::rk_f2f3(
    const grid_t &density_grid_aux_old, 
    grid_t &density_grid_aux_new, 
    grid_t &k_grid_out, 
    const double dt_in
)
{
    for (auto i=0; i<n_cells; i++)
    {
        k_grid_out[i] = nonlinear_rhs(i, density_grid_aux_old);
        density_grid_aux_new[i] = density_grid[i] + dt_in*k_grid_out[i];
    }
}

//! Runge-Kutta fourth function and stochastic step, all in the same loop
void BaseLangevin::rk_f4_and_stochastic(
    const grid_t &density_grid_aux_old, 
    const grid_t &k1_grid, 
    const grid_t &k2_grid, 
    const grid_t &k3_grid, 
    rng_t &rng
)
{
    mean_density = 0.0;
    for (auto i=0; i<n_cells; i++)
    {
        auto k4 = nonlinear_rhs(i, density_grid_aux_old);
        density_grid[i] += ( k1_grid[i] + 2*(k2_grid[i] + k3_grid[i]) + k4 )*(dt/6.0);
        #if !APPROXIMATE_POISSON_DISTBN
        poisson_rng = poisson_dist_t(lambda_product*density_grid[i]);
        gamma_rng = gamma_dist_t(poisson_rng(rng), 1.0);
        #else
        double mu = lambda_product * density_grid[i];
        if (mu > MU_THRESHOLD)
        {
            normal = normal_dist_t(mu, sqrt(mu));
            gamma_rng = gamma_dist_t(normal(rng), 1.0);
        }
        else
        {
            poisson_rng = poisson_dist_t(lambda_product*density_grid[i]);
            gamma_rng = gamma_dist_t(poisson_rng(rng), 1.0);
        }
        #endif
        density_grid[i] = gamma_rng(rng)/lambda;
        mean_density += density_grid[i];
    }    
    mean_density /= static_cast<double>(n_cells);    
}
