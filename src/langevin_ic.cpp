/**
 * @file langevin_ic.cpp
 * @brief Methods for setting up the initial condition of the Langevin model.
 */

#include "general_types.hpp"

bool BaseLangevin::initialize_grid(const Parameters p, rng_t& rng)
{
    int i_cell;
    switch (p.initial_condition)
    {
        case (InitialCondition::RANDOM_GAUSSIAN):
            ic_random_uniform(rng);
            return true;
        case (InitialCondition::CONSTANT_VALUE):
            ic_constant_value(p.ic_values.at(0));
            return true;
        case (InitialCondition::SINGLE_SEED):
            if (p.grid_dimension==GridDimension::D1)
            {
                i_cell = ( static_cast<int>(p.ic_values.at(1)) );
                if (i_cell<0 or i_cell>=p.n_x) { return false; }
            } 
            else if (p.grid_dimension==GridDimension::D2)
            {
                i_cell = (static_cast<int>(p.ic_values.at(1))
                        + static_cast<int>(p.ic_values.at(2))*p.n_x);
                if (i_cell<0 or i_cell>=p.n_x*p.n_y) { return false; }
            } 
            else if (p.grid_dimension==GridDimension::D3)
            {
                return false;
            } 
            else 
            { 
                return false; 
            }
            ic_single_seed(i_cell, p.ic_values.at(0));
            return true;
        case (InitialCondition::RANDOM_UNIFORM):
            ic_random_uniform(rng);
            return true;
        default:
            return false;
    }  
}

//! Set grid cells to have uniformly random values between min_value and max_value
void BaseLangevin::ic_random_uniform(
    rng_t &rng, 
    const double min_value, 
    const double max_value
)
{
    uniform_dist_t uniform(min_value, max_value);
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
    density_grid = grid_t(n_cells, density_value);
    mean_density = density_value;
}

//! Set all the grid cells to zero except a single specified cell
void BaseLangevin::ic_single_seed(const int i_cell, const double value)
{
    density_grid[i_cell] = value;
    mean_density = value / static_cast<double>(n_cells);
}    
