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

int SimDP::count_epochs() const
{
    // Count total number of time steps, just in case rounding causes problems
    int n_epochs;
    double t; 
    for (n_epochs=0, t=0; t<=p.t_final+p.dt; t+=p.dt, n_epochs++) {}
    return n_epochs;
}

bool SimDP::integrate(const int n_next_epochs)
{
    int i;
    double t; 
    // AAAARGHHH I don't yet know how to use a pointer to one
    //   of these integrate() functions; if I did, there would
    //   be no duplication here!
    // void (*integrate_fn)(rng_t&);
    void (DPLangevin::*integrate_fn)(rng_t&);
    integrate_fn = &DPLangevin::integrate_euler;
    // dpLangevin->*integrate_fn(*rng);
    // switch (p.integration_method)
    // {
    //     case (IntegrationMethod::EULER):
    //         integrate_fn = dpLangevin->integrate_euler(*rng);
    //         break;
    //     // case (IntegrationMethod::RUNGE_KUTTA):
    //     //     integrate_fn = dpLangevin->integrate_rungekutta;
    //     //     break;
    //     default:
    //         return false;
    // }

    if (epochs.size() < i_epoch+n_next_epochs)
    {
        std::cout << "Too many epochs: " 
            << epochs.size() << " < " << i_epoch+n_next_epochs << std::endl;
    }
    for (i=i_epoch, t=t_epoch; i<i_epoch+n_next_epochs; t+=p.dt, i++)
    {
        switch (p.integration_method)
        {
            case (IntegrationMethod::RUNGE_KUTTA): 
                dpLangevin->integrate_rungekutta(*rng);
                break;
            case (IntegrationMethod::EULER): 
                dpLangevin->integrate_euler(*rng);
                break;
            default:
                return false;
        }
        epochs[i] = t;
        mean_densities[i] = dpLangevin->get_mean_density();
    };
    i_epoch = i;
    t_epoch = t;
    return true;
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
    if (not (p.n_cells == p.n_x * p.n_y * p.n_z)) { 
        std::cout << "prep_density: grid size problem" << std::endl;
        return false; 
    }
    // Assume we're working with a 2d grid for now
    py_array_t density_array({p.n_x, p.n_y});
    auto density_proxy = density_array.mutable_unchecked();
    int i_x, i_y;
    for (auto i=0; i<p.n_cells; i++)
    {
        i_x = i % p.n_x;
        i_y = i / p.n_x;
        // std::cout << i << " " << i_x << " " << i_y << " " << std::endl;
        density_proxy(i_x, i_y) = dpLangevin->get_cell_density(i);
    };
    return_density = density_array;
    return true;
}