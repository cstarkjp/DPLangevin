/**
 * @file langevin_utils.cpp
 * @brief Utility methods to process the Langevin field grid.
 */

#include "general_core.hpp"

//! Make the cell density visible.
double Langevin::get_cell_density(const int i) const 
{
    return cell_density[i];
}

//! Make the density visible
double Langevin::get_mean_density() const
{
    // ! Make the density visible
    return mean_density;
}

//! Check the mean of the Poisson distribution
double Langevin::get_poisson_mean() const
{
    return lambda_product * mean_density;
}