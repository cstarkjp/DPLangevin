// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void Langevin::integrate_rungekutta(rng_t &rng)
{
    // Runge-Kutta integration of the non-linear term and diffusion.
    // Update of cells is done in the same loop as last RK step 
    // for efficiency
    rk_f1(aux_cell_old, k1);
    rk_f2f3(aux_cell_old, aux_cell_new, k2, dtm);
    // Swap contents is O(1), better than old = fast
    aux_cell_old.swap(aux_cell_new);            
    rk_f2f3(aux_cell_old, aux_cell_new, k3, dt);
    aux_cell_old.swap(aux_cell_new);
    rk_f4_and_stochastic(aux_cell_old, k1, k2, k3, rng);
}

// RK first function
void Langevin::rk_f1(dbl_vec_t &aux_cell, dbl_vec_t &k1)
{
    for (auto i=0; i<n_cells; i++)
    {
        k1[i] = nonlinear_rhs(i, cell_density);
        aux_cell[i] = cell_density[i] + dtm*k1[i];
    }
}

// RK second and third functions
void Langevin::rk_f2f3(
    const dbl_vec_t &aux_old, 
    dbl_vec_t &aux_new, 
    dbl_vec_t &k_out, 
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
void Langevin::rk_f4_and_stochastic(
    const dbl_vec_t &aux_old, 
    const dbl_vec_t &k1, 
    const dbl_vec_t &k2, 
    const dbl_vec_t &k3, 
    rng_t &rng
)
{
    mean_density = 0.0;
    for (auto i=0; i<n_cells; i++)
    {
        auto k4 = nonlinear_rhs(i, aux_old);
        cell_density[i] += dts*(k1[i] + 2*(k2[i] + k3[i]) + k4);
        #if !APPROXIMATE_POISSON_DISTBN
        poisson = int_poisson_dist_t(lambda_product*cell_density[i]);
        gamma = dbl_gamma_dist_t(poisson(rng), 1.0);
        #else
        double mu = lambda_product * cell_density[i];
        if (mu > MU_THRESHOLD)
        {
            normal = dbl_normal_dist_t(mu, sqrt(mu));
            gamma = dbl_gamma_dist_t(normal(rng), 1.0);
        }
        else
        {
            poisson = int_poisson_dist_t(lambda_product*cell_density[i]);
            gamma = dbl_gamma_dist_t(poisson(rng), 1.0);
        }
        #endif
        cell_density[i] = gamma(rng)/lambda;
        mean_density += cell_density[i];
    }    
    mean_density /= static_cast<double>(n_cells);    
}
