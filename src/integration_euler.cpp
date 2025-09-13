// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void LangevinBase::integrate_euler(RNG &rng)
{
    euler_and_stochastic(aux_cell_old, rng);
    cell_density.swap(aux_cell_old); 
}

void LangevinBase::euler_and_stochastic(dbl_vec_t &aux, RNG &rng)
{
    mean_density = 0.0;
    for (auto i=0; i<n_cells; i++)
    {
        double f = nonlinear_rhs(i, cell_density);
        aux[i] = cell_density[i] + dt*f;

        #if !APPROXIMATE_POISSON_DISTBN
        poisson = int_poisson_distbn(lambda_product * aux[i]);
        gamma = dbl_gamma_distbn(poisson(rng), 1.0);
        #else
        double mu = lambda_product * aux[i];
        if (mu > MU_THRESHOLD)
        {
            normal = dbl_normal_distbn(mu, sqrt(mu));
            gamma = dbl_gamma_distbn(normal(rng), 1.0);
        }
        else
        {
            poisson = int_poisson_distbn(lambda_product * aux[i]);
            gamma = dbl_gamma_distbn(poisson(rng), 1.0);
        }
        #endif
        aux[i]= gamma(rng)/lambda;
        mean_density += aux[i];
    }    
    mean_density /= static_cast<double>(n_cells);   
}