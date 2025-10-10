/**
 * @file application_dpsim_private.cpp
 * @brief Class to manage & run DPLangevin model simulation: private methods.
 */ 

#include "general_core.hpp"
#include "application_dpsim.hpp"

bool SimDP::construct_grid()
{
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
    int i_cell;
    switch (p.initial_condition)
    {
        case (InitialCondition::RANDOM_GAUSSIAN):
            dpLangevin->ic_random_uniform(*rng);
            return true;
        case (InitialCondition::CONSTANT_VALUE):
            dpLangevin->ic_constant_value(p.ic_values.at(0));
            return true;
        case (InitialCondition::SINGLE_SEED):
            if (p.grid_dimension==GridDimension::D1)
            {
                i_cell = ( static_cast<int>(p.ic_values.at(1)) );
                if (i_cell<0 or i_cell>=p.n_x) { return false; }
            } 
            else if (p.grid_dimension==GridDimension::D2)
            {
                i_cell = (static_cast<int>(p.ic_values.at(1))
                        + static_cast<int>(p.ic_values.at(2))*p.n_x);
                if (i_cell<0 or i_cell>=p.n_x*p.n_y) { return false; }
            } 
            else { return false; }
            dpLangevin->ic_single_seed(i_cell, p.ic_values.at(0));
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

bool SimDP::choose_integrator()
{
    switch (p.integration_method)
    {
        case (IntegrationMethod::RUNGE_KUTTA):
            integrator = &DPLangevin::integrate_rungekutta;
            return true;
        case (IntegrationMethod::EULER):
            integrator = &DPLangevin::integrate_euler;
            return true;
        default:
            return false;
    }
}

bool SimDP::integrate(const int n_next_epochs)
{
    // Check a further n_next_epochs won't exceed total permitted steps
    if (t_epochs.size() < i_next_epoch+n_next_epochs)
    {
        std::cout << "Too many epochs: " 
            << t_epochs.size() << " < " << i_next_epoch+n_next_epochs << std::endl;
        return false;
    }
    
    // Perform (possibly another another) n_next_epochs integration steps
    int i;
    double t; 
    // For the very first epoch, record mean density right now
    if (i_next_epoch==1) { 
        // dpLangevin->apply_boundary_conditions();
        mean_densities[0] = dpLangevin->get_mean_density(); 
        i_current_epoch = 0;
        t_current_epoch = 0;
    }
    // Loop over integration steps.
    // Effectively increment epoch counter and add to Δt to time counter
    // so that both point the state *after* each integration step is complete.
    // In so doing, we will record t_epochs.size() + 1 total integration steps.
    for (
        i=i_next_epoch, t=t_next_epoch; 
        i<i_next_epoch+n_next_epochs; 
        t+=p.dt, i++)
    {
        // Reapply boundary conditions prior to integrating
        dpLangevin->apply_boundary_conditions(p);
        // Perform a single integration over Δt
        (dpLangevin->*integrator)(*rng);
        // Record this epoch
        t_epochs[i] = t;
        mean_densities[i] = dpLangevin->get_mean_density();
        i_current_epoch = i;
        t_current_epoch = t;
    };
    // Set epoch and time counters to point to *after* the last integration step
    i_next_epoch = i;
    t_next_epoch = t;
    return true;
}

bool SimDP::prep_t_epochs()
{
    py_array_t epochs_array(n_epochs);
    auto epochs_proxy = epochs_array.mutable_unchecked();
    for (auto i=0; i<n_epochs; i++)
    {
        epochs_proxy(i) = t_epochs[i];
    };
    pyarray_t_epochs = epochs_array;
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
    pyarray_mean_densities = mean_densities_array;
    return true;
}

bool SimDP::prep_density_grid()
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
        density_proxy(i_x, i_y) = dpLangevin->get_density_grid_value(i);
    };
    pyarray_density = density_array;
    return true;
}