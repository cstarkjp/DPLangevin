// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

// Set cells to have uniformly random values between min_value and max_value
void DornicBase::random_intial_condition(
    RNG &rng, 
    const double min_value, 
    const double max_value
)
{
    dbl_uniform_distbn uniform(min_value, max_value);
    mean_density = 0.0;
    for (auto i=0; i<cell_density.size(); i++)
    {
        cell_density[i] += uniform(rng);
        mean_density += cell_density[i];
    }
    mean_density /= static_cast<double>(n_cells);
}

// Set all the cells to have same value
void DornicBase::homogeneous_initial_condition(const double density_value = 1.0)
{
    cell_density = dbl_vector(n_cells, density_value);
    mean_density = density_value;
}

// Set all the cells to zero except a single specified cell
void DornicBase::single_seed(const int i_node, const double value = 1.0)
{
    cell_density[i_node] = value;
    mean_density = value / static_cast<double>(n_cells);
}    
