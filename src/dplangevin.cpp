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
    n_cells = p.n_cells;
    density_grid = grid_t(n_cells, 0.0); 
    aux_grid2 = grid_t(n_cells);
    aux_grid1 = grid_t(n_cells);
    k1_grid = grid_t(n_cells, 0.0);
    k2_grid = grid_t(n_cells, 0.0);
    k3_grid = grid_t(n_cells, 0.0);
    dt = p.dt;
    dx = p.dx;
}

//! Method to set nonlinear coefficients in DP Langevin equation 
//! for deterministic integration step
void DPLangevin::set_nonlinear_coefficients(const Coefficients &coefficients)
{
    quadratic_coeff = coefficients.quadratic;
    D = coefficients.diffusion / (dx*dx);
}

//! Method to set nonlinear RHS of DP Langevin equation 
//! for deterministic integration step
double DPLangevin::nonlinear_rhs(const int i_cell, const grid_t &field) const
{
    // Non-linear term, which is quadratic in the DP Langevin
    const double quadratic_term 
        = -quadratic_coeff*field[i_cell]*field[i_cell];

    // Integration of diffusion
    double diffusion_sum = 0.0;
    auto n_neighbor_cells = grid_wiring[i_cell].size();
    for (auto i=0; i<n_neighbor_cells; i++)
    {
        auto i_neighbor = grid_wiring[i_cell][i];
        diffusion_sum += field[i_neighbor];
    }
    diffusion_sum = D*(diffusion_sum - n_neighbor_cells*field[i_cell]);
    return diffusion_sum + quadratic_term;
}
