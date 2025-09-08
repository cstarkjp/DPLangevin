// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include <pybind11/numpy.h>
#include "base.hpp"
#include "applications.hpp"

auto dp_demo(
    double linear, double quadratic, double diffusion, double noise, 
    int n_cells, double t_max, double dx, double dt, int random_seed
) -> results_t
{
    // Set up parameters, coefficients etc
    Coefficients f_coeffs (linear, quadratic, diffusion, noise);
    Parameters parameters (n_cells, t_max, dx, dt, random_seed);
    RNG rng(random_seed); 
    f_coeffs.print();
    parameters.print();

    // Initialize
    Dornic dornic(parameters);
    dornic.construct_2D_grid();
    dornic.set_coefficients(f_coeffs);
    dornic.random_intial_condition(rng);

    // Integrate
    int n_epochs, i;
    double t; 
    // Count total number of time steps, just in case rounding causes problems
    for (n_epochs=0, t=0; t<=t_max+dt; t+=dt, n_epochs++) {}
    std::cout << "Total number of epochs: " << n_epochs << std::endl;
    dbl_vector epochs(n_epochs, 0.0);
    dbl_vector mean_densities(n_epochs, 0.0);
    for (i=0, t=0; i<n_epochs; t+=dt, i++)
    {
        dornic.integration(rng);
        epochs[i] = t;
        mean_densities[i] = dornic.density();
    }
    
    // Prepare return array
    results_t results({n_epochs, 2});
    auto ra = results.mutable_unchecked();
    for (auto i=0; i<n_epochs; i++)
    {
        ra(i, 0) = epochs[i];
        ra(i, 1) = mean_densities[i];
    };
    return results;
} 