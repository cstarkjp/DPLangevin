/**
 * @file langevin_integrate_rungekutta.cpp
 * @brief Methods to carry out 4th-order Runge-Kutta integration.
 */

#include "langevin_types.hpp"
#include "langevin_base.hpp"

//! Runge-Kutta integration of the nonlinear and diffusion terms 
//! in the Langevin equation.
//! Update of cells is done in the same loop as last Runge-Kutta step 
//! for efficiency.
void BaseLangevin::integrate_rungekutta(rng_t& rng)
{
    auto rk_f1 = [&](
        grid_t& aux_grid, grid_t& k1_grid, 
        const double dt_fraction
    )
    {
        for (auto i=0; i<n_cells; i++)
        {
            k1_grid[i] = nonlinear_rhs(i, density_grid);
            aux_grid[i] = density_grid[i] + dt_fraction*k1_grid[i];
        }
    };
    auto rk_f2f3 = [&](
        const grid_t& aux_grid_in, grid_t& aux_grid_out, grid_t& k23_grid, 
        const double dt_fraction
    )
    {
        for (auto i=0; i<n_cells; i++)
        {
            k23_grid[i] = nonlinear_rhs(i, aux_grid_in);
            aux_grid_out[i] = density_grid[i] + dt_fraction*k23_grid[i];
        }
    };
    auto rk_f4_stochastic = [&](
        const grid_t& aux_grid, const grid_t& k1_grid, const grid_t& k2_grid, 
        const grid_t& k3_grid, rng_t &rng, const double dt_fraction
    )
    {
        mean_density = 0.0;
        for (auto i=0; i<n_cells; i++)
        {
            // Runge-Kutta 4th step
            auto k4 = nonlinear_rhs(i, aux_grid);
            density_grid[i] 
                += ( k1_grid[i] + 2*(k2_grid[i] + k3_grid[i]) + k4 )*dt_fraction;
            // Stochastic step
            poisson_rng = poisson_dist_t(lambda_product*density_grid[i]);
            gamma_rng = gamma_dist_t(poisson_rng(rng), 1.0);
            density_grid[i] = gamma_rng(rng)/lambda;
            // Incrementally compute mean density
            mean_density += density_grid[i];
        }    
        mean_density /= static_cast<double>(n_cells);    
    };

    rk_f1(aux_grid1, k1_grid, dt/2);
    rk_f2f3(aux_grid1, aux_grid2, k2_grid, dt/2);
    rk_f2f3(aux_grid2, aux_grid1, k3_grid, dt);
    rk_f4_stochastic(aux_grid1, k1_grid, k2_grid, k3_grid, rng, dt/6);
}

// //! Runge-Kutta first step
// void BaseLangevin::rk_f1(
//     grid_t& aux_grid, 
//     grid_t& k1_grid, 
//     const double dt_fraction
// )
// {
//     for (auto i=0; i<n_cells; i++)
//     {
//         k1_grid[i] = nonlinear_rhs(i, density_grid);
//         aux_grid[i] = density_grid[i] + dt_fraction*k1_grid[i];
//     }
// }

// //! Runge-Kutta second and third steps
// void BaseLangevin::rk_f2f3(
//     const grid_t& aux_grid_in, 
//     grid_t& aux_grid_out, 
//     grid_t& k23_grid, 
//     const double dt_fraction
// )
// {
//     for (auto i=0; i<n_cells; i++)
//     {
//         k23_grid[i] = nonlinear_rhs(i, aux_grid_in);
//         aux_grid_out[i] = density_grid[i] + dt_fraction*k23_grid[i];
//     }
// }

// //! Runge-Kutta fourth step and stochastic step, all at once
// void BaseLangevin::rk_f4_stochastic(
//     const grid_t& aux_grid, 
//     const grid_t& k1_grid, 
//     const grid_t& k2_grid, 
//     const grid_t& k3_grid, 
//     rng_t &rng,
//     const double dt_fraction
// )
// {
//     mean_density = 0.0;
//     for (auto i=0; i<n_cells; i++)
//     {
//         // Runge-Kutta 4th step
//         auto k4 = nonlinear_rhs(i, aux_grid);
//         density_grid[i] 
//             += ( k1_grid[i] + 2*(k2_grid[i] + k3_grid[i]) + k4 )*dt_fraction;
//         // Stochastic step
//         poisson_rng = poisson_dist_t(lambda_product*density_grid[i]);
//         gamma_rng = gamma_dist_t(poisson_rng(rng), 1.0);
//         density_grid[i] = gamma_rng(rng)/lambda;
//         // Incrementally compute mean density
//         mean_density += density_grid[i];
//     }    
//     mean_density /= static_cast<double>(n_cells);    
// }
