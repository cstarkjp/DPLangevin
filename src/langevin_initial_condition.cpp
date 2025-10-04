/**
 * @file langevin_initial_condition.cpp
 * @brief Methods for setting up the initial condition of the Langevin model.
 */

#include "general_core.hpp"

//! Set grid cells to have uniformly random values between min_value and max_value
void BaseLangevin::ic_random_uniform(
    rng_t &rng, 
    const double min_value, 
    const double max_value
)
{
    dbl_uniform_dist_t uniform(min_value, max_value);
    mean_density = 0.0;
    for (auto i=0; i<density_grid.size(); i++)
    {
        density_grid[i] = uniform(rng);
        mean_density += density_grid[i];
    }
    mean_density /= static_cast<double>(n_cells);
}

//! Set all the grid cells to have same value
void BaseLangevin::ic_constant_value(const double density_value)
{
    density_grid = dbl_vec_t(n_cells, density_value);
    mean_density = density_value;
}

//! Set all the grid cells to zero except a single specified cell
void BaseLangevin::ic_single_seed(const int i_node, const double value)
{
    density_grid[i_node] = value;
    mean_density = value / static_cast<double>(n_cells);
}    
