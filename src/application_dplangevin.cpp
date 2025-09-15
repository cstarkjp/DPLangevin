// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include <pybind11/numpy.h>
#include <string>
#include "general_core.hpp"
#include "application_dplangevin.hpp"

DPLangevin::DPLangevin(Parameters p)
{
    n_cells = p.n_cells;
    cell_density = dbl_vec_t(n_cells, 0.0); 
    aux_cell_new = dbl_vec_t(n_cells);
    aux_cell_old = dbl_vec_t(n_cells);
    k1 = dbl_vec_t(n_cells, 0.0);
    k2 = dbl_vec_t(n_cells, 0.0);
    k3 = dbl_vec_t(n_cells, 0.0);
    k4 = dbl_vec_t(n_cells, 0.0);
    dt = p.dt;
    dx = p.dx;
    dtm = 0.5*dt;
    dts = dt/6.0;
}

void DPLangevin::set_nonlinear_coefficients(const Coefficients &f_coefficients)
{
    quadratic_coeff = f_coefficients.quadratic;
    D = f_coefficients.diffusion / (dx*dx);
}

double DPLangevin::nonlinear_rhs(const int i_cell, const dbl_vec_t &field) const
{
    // Non-linear terms
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
