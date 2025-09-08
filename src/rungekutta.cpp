// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

// RK first function
void DornicBase::f1(dbl_vector &aux_cell, dbl_vector &k1)
{
    for (auto i=0; i<n_cells; i++)
    {
        k1[i] = nonlinear_rhs(i, cell_density);
        aux_cell[i] = cell_density[i] + dtm*k1[i];
    }
}

// RK second and third functions
void DornicBase::f2f3(
    const dbl_vector &aux_old, 
    dbl_vector &aux_new, 
    dbl_vector &k_out, 
    const double dt_in
)
{
    for (auto i=0; i<n_cells; i++)
    {
        k_out[i] = nonlinear_rhs(i, aux_old);
        aux_new[i] = cell_density[i] + dt_in*k_out[i];
    }
}

// RK fourth function and stochastic step, all in the same loop
void DornicBase::f4_and_stochastic(
    const dbl_vector &aux_old, 
    const dbl_vector &k1, 
    const dbl_vector &k2, 
    const dbl_vector &k3, 
    RNG &rng
)
{
    mean_density = 0.0;
    for (auto i=0; i<n_cells; i++)
    {
        auto k4 = nonlinear_rhs(i, aux_old);
        cell_density[i] += dts*(k1[i] + 2*(k2[i] + k3[i]) + k4);
        #if !APPROXIMATE_POISSON_DISTBN
        poisson = int_poisson_distbn(lambda_product*cell_density[i]);
        gamma = dbl_gamma_distbn(poisson(rng), 1.0);
        #else
        double mu = lambda_product * cell_density[i];
        if (mu > MU_THRESHOLD)
        {
            normal = dbl_normal_distbn(mu, sqrt(mu));
            gamma = dbl_gamma_distbn(normal(rng), 1.0);
        }
        else
        {
            poisson = int_poisson_distbn(lambda_product*cell_density[i]);
            gamma = dbl_gamma_distbn(poisson(rng), 1.0);
        }
        #endif
        cell_density[i] = gamma(rng)/lambda;
        mean_density += cell_density[i];
    }    
    mean_density /= static_cast<double>(n_cells);    
}
