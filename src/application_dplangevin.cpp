/**
 * @file application_dplangevin.cpp
 * @brief Redefinition of BaseLangevin constructor; implementation of stub methods.
 */

#include <pybind11/numpy.h>
#include <string>
#include "general_core.hpp"
#include "application_dplangevin.hpp"

/**
 * @brief Redefinition of BaseLangevin class constructor
 */
DPLangevin::DPLangevin(Parameters p)
{
    n_cells = p.n_cells;
    density_grid = dbl_vec_t(n_cells, 0.0); 
    density_grid_aux_new = dbl_vec_t(n_cells);
    density_grid_aux_old = dbl_vec_t(n_cells);
    k1 = dbl_vec_t(n_cells, 0.0);
    k2 = dbl_vec_t(n_cells, 0.0);
    k3 = dbl_vec_t(n_cells, 0.0);
    k4 = dbl_vec_t(n_cells, 0.0);
    dt = p.dt;
    dx = p.dx;
    dtm = 0.5*dt;
    dts = dt/6.0;
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
double DPLangevin::nonlinear_rhs(const int i_cell, const dbl_vec_t &field) const
{
    // Non-linear term, which is quadratic in the DP Langevin
    const double quadratic_term 
        = -quadratic_coeff*field[i_cell]*field[i_cell];

    // Integration of diffusion
    double diffusion_sum = 0.0;
    int n_neighbors = neighbors[i_cell].size();
    for (auto i=0; i<n_neighbors; i++)
    {
        auto i_neighbor = neighbors[i_cell][i];
        diffusion_sum += field[i_neighbor];
    }
    diffusion_sum = D*(diffusion_sum - n_neighbors*field[i_cell]);
    return diffusion_sum + quadratic_term;
}
