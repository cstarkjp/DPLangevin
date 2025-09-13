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

    auto nonlinear_rhs(const int i_cell, const dbl_vec_t &field) 
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
    dbl_vec_t& epochs, dbl_vec_t& mean_densities
)
{
    int i;
    double t; 

    switch (parameters.integration_method)
    {
        case (IntegrationMethod::EULER):
            for (i=0, t=0; i<epochs.size(); t+=parameters.dt, i++)
            {
                dpLangevin.integrate_euler(rng);
                epochs[i] = t;
                mean_densities[i] = dpLangevin.get_mean_density();
            };
            break;
        case (IntegrationMethod::RUNGE_KUTTA):
        default:
            for (i=0, t=0; i<epochs.size(); t+=parameters.dt, i++)
            {
                dpLangevin.integrate_rungekutta(rng);
                epochs[i] = t;
                mean_densities[i] = dpLangevin.get_mean_density();
            };
            break;
    }
}

results_t prepare_return_array(
    const dbl_vec_t& cell_density,
    const dbl_vec_t& epochs, const dbl_vec_t& mean_densities
)
{
    results_t results({static_cast<int>(epochs.size()), 2});
    auto array_proxy = results.mutable_unchecked();
    for (auto i=0; i<epochs.size(); i++)
    {
        array_proxy(i, 0) = epochs[i];
        array_proxy(i, 1) = mean_densities[i];
    };
    return results;
}


auto dp(
    const double linear, const double quadratic, 
    const double diffusion, const double noise, 
    const double t_max, const double dx, const double dt, const int random_seed,
    const GridDimension grid_dimension, 
    const int_vec_t& grid_size,
    const GridTopology grid_topology,
    const BoundaryCondition boundary_condition,
    const InitialCondition initial_condition,
    const IntegrationMethod integration_method
) -> Results
{
    Coefficients f_coeffs (linear, quadratic, diffusion, noise);
    Parameters parameters (
        t_max, dx, dt, random_seed,
        grid_dimension, grid_size, grid_topology, 
        boundary_condition, initial_condition, integration_method
    );
    RNG rng(parameters.random_seed); 
    DPLangevin dpLangevin(parameters);
    f_coeffs.print();
    parameters.print();

    construct_grid(dpLangevin, parameters);
    initialize_grid(dpLangevin, parameters, rng);
    dpLangevin.set_coefficients(f_coeffs);
    auto n_epochs = count_epochs(parameters);
    dbl_vec_t epochs(n_epochs, 0.0);
    dbl_vec_t mean_densities(n_epochs, 0.0);
    // std::array<double, n_epochs> epochs;
    // std::array<double, n_epochs> mean_densities;
    integrate(
        dpLangevin, parameters, rng, epochs, mean_densities
    );
    
    // return prepare_return_array(
    //     dpLangevin.get_density(), epochs, mean_densities
    // );
    Results results(
        n_epochs, 
        parameters.n_cells, 
        parameters.n_x, parameters.n_y, parameters.n_z
    );
    results.prep_epochs(epochs);
    results.prep_mean_densities(mean_densities);
    results.prep_density(dpLangevin.get_density());
    return results;
} 