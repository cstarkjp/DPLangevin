// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "general_core.hpp"
#include "application_dpsim.hpp"

bool SimDP::construct_grid()
{
    // std::cout << "construct_grid::  dpLangevin = " << dpLangevin << std::endl;
    switch (p.grid_dimension)
    {
        case (GridDimension::D1):
            return dpLangevin->construct_1D_grid(p);
        case (GridDimension::D2):
            return dpLangevin->construct_2D_grid(p);
        default:
            return false;
    }    
}

bool SimDP::initialize_grid()
    {
        switch (p.initial_condition)
        {
            case (InitialCondition::RANDOM_GAUSSIAN):
                dpLangevin->ic_random_uniform(*rng);
                return true;
            case (InitialCondition::CONSTANT_VALUE):
                dpLangevin->ic_constant_value(1.0);
                return true;
            case (InitialCondition::SINGLE_SEED):
                // n_cells is not yet set
                dpLangevin->ic_single_seed(
                    static_cast<double>(p.n_cells)/2.0, 1.0
                );
                return true;
            case (InitialCondition::RANDOM_UNIFORM):
                dpLangevin->ic_random_uniform(*rng);
                return true;
            default:
                return false;
        }  
    }

int SimDP::count_epochs()
{
    // Count total number of time steps, just in case rounding causes problems
    int n_epochs;
    double t; 
    for (n_epochs=0, t=0; t<=p.t_max+p.dt; t+=p.dt, n_epochs++) {}
    return n_epochs;
}

bool SimDP::integrate(dbl_vec_t& epochs, dbl_vec_t& mean_densities)
{
    int i;
    double t; 

    // std::cout << "integrate::  n_epochs = " << n_epochs << std::endl;
    // std::cout << "integrate::  epochs.size = " << epochs.size() << std::endl;
    // std::cout << "integrate::  mean_densities.size = " << mean_densities.size() << std::endl;
    switch (p.integration_method)
    {
        case (IntegrationMethod::EULER):
            // std::cout << "integrate::  Euler "<< std::endl;
            for (i=0, t=0; i<epochs.size(); t+=p.dt, i++)
            {
                dpLangevin->integrate_euler(*rng);
                epochs[i] = t;
                mean_densities[i] = dpLangevin->get_mean_density();
            };
            return true;
        case (IntegrationMethod::RUNGE_KUTTA):
            // std::cout << "integrate::  Runge-Kutta "<< std::endl;
            // std::cout << "integrate::  dpLangevin = " << dpLangevin << std::endl;
            for (i=0, t=0; i<epochs.size(); t+=p.dt, i++)
            {
                dpLangevin->integrate_rungekutta(*rng);
                epochs[i] = t;
                mean_densities[i] = dpLangevin->get_mean_density();
            };
            return true;
        default:
            return false;
    }
}

bool SimDP::prep_epochs()
{
    py_array_t epochs_array(n_epochs);
    auto epochs_proxy = epochs_array.mutable_unchecked();
    for (auto i=0; i<n_epochs; i++)
    {
        epochs_proxy(i) = epochs[i];
    };
    return_epochs = epochs_array;
    return true;
}

bool SimDP::prep_mean_densities()
{
    py_array_t mean_densities_array(n_epochs);
    auto mean_densities_proxy = mean_densities_array.mutable_unchecked();
    for (auto i=0; i<n_epochs; i++)
    {
        mean_densities_proxy(i) = mean_densities[i];
    };
    return_mean_densities = mean_densities_array;
    return true;
}

bool SimDP::prep_density()
{
    py_array_t density_array({p.n_x, p.n_y});
    // auto density_proxy = density_array.mutable_unchecked();
    // std::cout << "prep_density: p.n_cells = " << p.n_cells << std::endl;
    // std::cout << "prep_density: p.n_x = " << p.n_x << std::endl;
    // std::cout << "prep_density: p.n_y = " << p.n_y << std::endl;
    // std::cout << "prep_density: p.n_z = " << p.n_z << std::endl;
    if (not (p.n_cells == p.n_x * p.n_y * p.n_z)) { 
        std::cout << "prep_density: failed" << std::endl;
        return false; 
    }
    for (auto i=0; i<p.n_cells; i++)
    {
        // density_proxy(i, 0) = epochs[i];
        // density_proxy(i, 1) = mean_densities[i];
    };
    return_density = density_array;
    return true;
}