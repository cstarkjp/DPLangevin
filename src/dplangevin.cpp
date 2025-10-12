/**
 * @file dplangevin.cpp
 * @brief Redefinition of BaseLangevin constructor; implementation of stub methods.
 */

#include <pybind11/numpy.h>
#include <string>
#include "dplangevin.hpp"

/**
 * @brief Redefinition of BaseLangevin class constructor
 */
DPLangevin::DPLangevin(Parameters p)
{
    // "Local" copies
    n_cells = p.n_cells;
    dt = p.dt;
    dx = p.dx;
    // The Langevin density field grid as a 1d vector
    density_grid = grid_t(n_cells, 0.0); 
    // Supplementary grids (as 1d vectors of same length)
    aux_grid1 = grid_t(n_cells, 0.0);
    aux_grid2 = grid_t(n_cells, 0.0);
    k1_grid = grid_t(n_cells, 0.0);
    k2_grid = grid_t(n_cells, 0.0);
    k3_grid = grid_t(n_cells, 0.0);
}

//! Method to set nonlinear coefficients in DP Langevin equation 
//! for deterministic integration step
void DPLangevin::set_nonlinear_coefficients(const Coefficients& coefficients)
{
    quadratic_coefficient = coefficients.quadratic;
    diffusion_coefficient = coefficients.diffusion / (dx*dx);
}

//! Method to set nonlinear RHS of DP Langevin equation 
//! for deterministic integration step
double DPLangevin::nonlinear_rhs(const int i_cell, const grid_t& density) const
{
    // Non-linear term, which is quadratic in the DP Langevin
    const double quadratic_term 
        = -quadratic_coefficient*density[i_cell]*density[i_cell];

    // Integration of diffusion
    double diffusion_sum = 0.0;
    auto cell_wiring = grid_wiring[i_cell];
    auto n_neighbors = cell_wiring.size();
    for (auto j_wire=0; j_wire<n_neighbors; j_wire++)
    {
        auto j_neighbor_cell = cell_wiring[j_wire];
        diffusion_sum += density[j_neighbor_cell];
    }
    auto diffusion_term = (
        diffusion_coefficient*(diffusion_sum - n_neighbors*density[i_cell])
    );
    return diffusion_term + quadratic_term;
}
