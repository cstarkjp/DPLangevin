// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include <pybind11/numpy.h>
#include "base.hpp"
#include "applications.hpp"

class Dornic_DP : public Dornic {

public:
    double quadratic_coeff;
    double D;

    Dornic_DP(Parameters params) : Dornic(params) {}

    void set_nonlinear_coefficients(const Coefficients &f_coefficients)
    override
    {
        quadratic_coeff = f_coefficients.quadratic;
        D = f_coefficients.diffusion / (dx*dx);
    }

    double nonlinear_rhs(const int i_node, const dbl_vector &field) 
    const
    override
    {
        // Non-linear terms
        const double quadratic_term = -quadratic_coeff*field[i_node]*field[i_node];

        // Integration of diffusion
        double diffusion_sum = 0.0;
        int n_neighbors = neighbors[i_node].size();
        for (auto i=0; i<n_neighbors; i++)
        {
            auto i_neighbor = neighbors[i_node][i];
            diffusion_sum += field[i_neighbor];
        }
        diffusion_sum = D*(diffusion_sum - n_neighbors*field[i_node]);
        return diffusion_sum + quadratic_term;
    }

    void check(void) const
    {
        std::cout << "In Dornic_DP class instance" << std::endl;
    }
};

auto dp(
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
    Dornic_DP dornic(parameters);
    dornic.check();
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