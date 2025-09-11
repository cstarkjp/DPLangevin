// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void DornicBase::integrate(RNG &rng)
{
    #if (USE_RUNGEKUTTA)
    // Runge-Kutta integration of the non-linear term and diffusion.
    // Update of cells is done in the same loop as last RK step 
    // for efficiency
    f1(aux_cell_old, k1);
    f2f3(aux_cell_old, aux_cell_new, k2, dtm);
    // Swap contents is O(1), better than old = fast
    aux_cell_old.swap(aux_cell_new);            
    f2f3(aux_cell_old, aux_cell_new, k3, dt);
    aux_cell_old.swap(aux_cell_new);
    f4_and_stochastic(aux_cell_old, k1, k2, k3, rng);
    #else
    euler_and_stochastic(aux_cell_old, rng);
    cell_density.swap(aux_cell_old); 
    #endif
}

void DornicBase::integrate_rungekutta(RNG &rng)
{
    // Runge-Kutta integration of the non-linear term and diffusion.
    // Update of cells is done in the same loop as last RK step 
    // for efficiency
    f1(aux_cell_old, k1);
    f2f3(aux_cell_old, aux_cell_new, k2, dtm);
    // Swap contents is O(1), better than old = fast
    aux_cell_old.swap(aux_cell_new);            
    f2f3(aux_cell_old, aux_cell_new, k3, dt);
    aux_cell_old.swap(aux_cell_new);
    f4_and_stochastic(aux_cell_old, k1, k2, k3, rng);
}

void DornicBase::integrate_euler(RNG &rng)
{
    euler_and_stochastic(aux_cell_old, rng);
    cell_density.swap(aux_cell_old); 
}


