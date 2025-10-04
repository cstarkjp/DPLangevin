/**
 * @file langevin_integration_rungekutta.cpp
 * @brief Methods to carry out 4th-order Runge-Kutta integration.
 */

#include "general_core.hpp"

//! Runge-Kutta integration of the nonlinear and diffusion terms 
//! in the Langevin equation.
//! Update of cells is done in the same loop as last Runge-Kutta step 
//! for efficiency.
void BaseLangevin::integrate_rungekutta(rng_t &rng)
{

    rk_f1(density_grid_aux_old, k1);
    rk_f2f3(density_grid_aux_old, density_grid_aux_new, k2, dtm);
    // Swap contents is O(1), better than old = fast
    density_grid_aux_old.swap(density_grid_aux_new);            
    rk_f2f3(density_grid_aux_old, density_grid_aux_new, k3, dt);
    density_grid_aux_old.swap(density_grid_aux_new);
    rk_f4_and_stochastic(density_grid_aux_old, k1, k2, k3, rng);
}

//! Runge-Kutta first function
void BaseLangevin::rk_f1(dbl_vec_t &density_grid_aux, dbl_vec_t &k1)
{
    for (auto i=0; i<n_cells; i++)
    {
        k1[i] = nonlinear_rhs(i, density_grid);
        density_grid_aux[i] = density_grid[i] + dtm*k1[i];
    }
}

//! Runge-Kutta second and third functions
void BaseLangevin::rk_f2f3(
    const dbl_vec_t &density_grid_aux_old, 
    dbl_vec_t &density_grid_aux_new, 
    dbl_vec_t &k_out, 
    const double dt_in
)
{
    for (auto i=0; i<n_cells; i++)
    {
        k_out[i] = nonlinear_rhs(i, density_grid_aux_old);
        density_grid_aux_new[i] = density_grid[i] + dt_in*k_out[i];
    }
}

//! Runge-Kutta fourth function and stochastic step, all in the same loop
void BaseLangevin::rk_f4_and_stochastic(
    const dbl_vec_t &density_grid_aux_old, 
    const dbl_vec_t &k1, 
    const dbl_vec_t &k2, 
    const dbl_vec_t &k3, 
    rng_t &rng
)
{
    mean_density = 0.0;
    for (auto i=0; i<n_cells; i++)
    {
        auto k4 = nonlinear_rhs(i, density_grid_aux_old);
        density_grid[i] += dts*(k1[i] + 2*(k2[i] + k3[i]) + k4);
        #if !APPROXIMATE_POISSON_DISTBN
        poisson_rng = int_poisson_dist_t(lambda_product*density_grid[i]);
        gamma_rng = dbl_gamma_dist_t(poisson_rng(rng), 1.0);
        #else
        double mu = lambda_product * density_grid[i];
        if (mu > MU_THRESHOLD)
        {
            normal = dbl_normal_dist_t(mu, sqrt(mu));
            gamma_rng = dbl_gamma_dist_t(normal(rng), 1.0);
        }
        else
        {
            poisson_rng = int_poisson_dist_t(lambda_product*density_grid[i]);
            gamma_rng = dbl_gamma_dist_t(poisson_rng(rng), 1.0);
        }
        #endif
        density_grid[i] = gamma_rng(rng)/lambda;
        mean_density += density_grid[i];
    }    
    mean_density /= static_cast<double>(n_cells);    
}
