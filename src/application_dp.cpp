// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include <pybind11/numpy.h>
#include <string>
#include "core.hpp"
#include "application_dp.hpp"

class DPLangevin : public LangevinBase 
{
public:
    double quadratic_coeff;
    double D;

    DPLangevin(Parameters params) : LangevinBase(params) {}

    void set_nonlinear_coefficients(const Coefficients &f_coefficients)
    override
    {
        quadratic_coeff = f_coefficients.quadratic;
        D = f_coefficients.diffusion / (dx*dx);
    }

    auto nonlinear_rhs(const int i_cell, const dbl_vector &field) 
    const -> double
    override
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
};

void construct_grid(DPLangevin& dpLangevin, const Parameters parameters)
{
    switch (parameters.grid_dimension)
    {
        case (GridDimension::D1):
            dpLangevin.construct_1D_grid(parameters);
            break;
        case (GridDimension::D2):
        default:
            dpLangevin.construct_2D_grid(parameters);
            break;
    }    
}

void initialize_grid(DPLangevin& dpLangevin, const Parameters parameters, RNG &rng)
{
    switch (parameters.initial_condition)
    {
        case (InitialCondition::RANDOM_GAUSSIAN):
            dpLangevin.ic_random_uniform(rng);
            break;
        case (InitialCondition::CONSTANT_VALUE):
            dpLangevin.ic_constant_value(1.0);
            break;
        case (InitialCondition::SINGLE_SEED):
            dpLangevin.ic_single_seed(parameters.n_cells/2, 1.0);
            break;
        case (InitialCondition::RANDOM_UNIFORM):
        default:
            dpLangevin.ic_random_uniform(rng);
            break;
    }  
}

int count_epochs(const Parameters parameters)
{
    int n_epochs;
    double t; 
    // Count total number of time steps, just in case rounding causes problems
    for (
        n_epochs=0, t=0; 
        t<=parameters.t_max+parameters.dt; 
        t+=parameters.dt, n_epochs++
    ) {}
    return n_epochs;
}

void integrate(
    DPLangevin& dpLangevin, const Parameters parameters, RNG& rng,
    const int n_epochs, dbl_vector& epochs, dbl_vector& mean_densities
)
{
    int i;
    double t; 

    switch (parameters.integration_method)
    {
        case (IntegrationMethod::EULER):
            for (i=0, t=0; i<n_epochs; t+=parameters.dt, i++)
            {
                dpLangevin.integrate_euler(rng);
                epochs[i] = t;
                mean_densities[i] = dpLangevin.density();
            };
            break;
        case (IntegrationMethod::RUNGE_KUTTA):
        default:
            for (i=0, t=0; i<n_epochs; t+=parameters.dt, i++)
            {
                dpLangevin.integrate_rungekutta(rng);
                epochs[i] = t;
                mean_densities[i] = dpLangevin.density();
            };
            break;
    }
}

results_t prepare_return_array(
    const int n_epochs, 
    const dbl_vector& epochs, 
    const dbl_vector& mean_densities
)
{
    results_t results({n_epochs, 2});
    auto ra = results.mutable_unchecked();
    for (auto i=0; i<n_epochs; i++)
    {
        ra(i, 0) = epochs[i];
        ra(i, 1) = mean_densities[i];
    };
    return results;
}

results_t dp(
    const double linear, const double quadratic, 
    const double diffusion, const double noise, 
    const double t_max, const double dx, const double dt, const int random_seed,
    const GridDimension grid_dimension, 
    const int_vector& grid_size,
    const GridTopology grid_topology,
    const BoundaryCondition boundary_condition,
    const InitialCondition initial_condition,
    const IntegrationMethod integration_method
)
{
    Coefficients f_coeffs (linear, quadratic, diffusion, noise);
    Parameters parameters (
        t_max, dx, dt, random_seed,
        grid_dimension, 
        grid_size, 
        grid_topology, 
        boundary_condition, 
        initial_condition, 
        integration_method
    );
    RNG rng(parameters.random_seed); 
    f_coeffs.print();
    parameters.print();

    DPLangevin dpLangevin(parameters);
    construct_grid(dpLangevin, parameters);
    initialize_grid(dpLangevin, parameters, rng);
    dpLangevin.set_coefficients(f_coeffs);
    int n_epochs = count_epochs(parameters);
    dbl_vector epochs(n_epochs, 0.0), mean_densities(n_epochs, 0.0);
    integrate(
        dpLangevin, parameters, rng, n_epochs, epochs, mean_densities
    );
    
    return prepare_return_array(n_epochs, epochs, mean_densities);
} 