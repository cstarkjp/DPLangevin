// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "base.hpp"

// // User-defined variables
// double quadratic_coeff;
// double D;

void Dornic::set_nonlinear_coefficients(const Coefficients &f_coefficients)
{
    // quadratic_coeff = f_coefficients.quadratic;
    // D = f_coefficients.diffusion / (dx*dx);
}

auto Dornic::nonlinear_rhs(const int i_node, const dbl_vector &field) 
    const -> double
{
    // // Non-linear terms
    // const double quadratic_term = -quadratic_coeff*field[i_node]*field[i_node];

    // // Integration of diffusion
    // double diffusion_sum = 0.0;
    // int n_neighbors = neighbors[i_node].size();
    // for (auto i=0; i<n_neighbors; i++)
    // {
    //     auto i_neighbor = neighbors[i_node][i];
    //     diffusion_sum += field[i_neighbor];
    // }
    // diffusion_sum = D*(diffusion_sum - n_neighbors*field[i_node]);
    // return diffusion_sum + quadratic_term;
    return 0.0;
}
