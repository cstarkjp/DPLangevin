/**
 * @file application_dpsim_public.cpp
 * @brief Class to manage & run DPLangevin model simulation: public methods.
 */ 

#include "general_core.hpp"
#include "application_dpsim.hpp"

SimDP::SimDP(
    const double linear, const double quadratic,
    const double diffusion, const double noise, 
    const double t_final, 
    const double dx, const double dt, 
    const int random_seed,
    const GridDimension grid_dimension,
    const int_vec_t grid_size,
    const gt_vec_t grid_topologies,
    const bc_vec_t boundary_conditions,
    const dbl_vec_t bc_values,
    const InitialCondition initial_condition,
    const dbl_vec_t ic_values,
    const IntegrationMethod integration_method
) : coefficients(linear, quadratic, diffusion, noise),
    p(
        t_final, 
        dx, 
        dt, 
        random_seed,
        grid_dimension, 
        grid_size, 
        grid_topologies,
        boundary_conditions,
        bc_values,
        initial_condition, 
        ic_values, 
        integration_method
    )
{
    rng = new rng_t(p.random_seed); 
    dpLangevin = new DPLangevin(p);
    coefficients.print();
    p.print();
}

bool SimDP::initialize()
{
    if (not construct_grid()) { 
        std::cout 
            << "SimDP::initialize failure: couldn't construct grid" 
            << std::endl;
        return false; 
    }
    if (not initialize_grid()) { 
        std::cout 
            << "SimDP::initialize failure: couldn't initialize grid" 
            << std::endl;
        return false; 
    }
    dpLangevin->set_coefficients(coefficients);
    n_epochs = count_epochs();
    t_epochs = dbl_vec_t(n_epochs, 0.0);
    mean_densities = dbl_vec_t(n_epochs, 0.0);
    // Treat epoch#0 as the initial grid state
    // So after initialization, we are nominally at epoch#1
    i_next_epoch = 1;
    t_next_epoch = p.dt;
    if (not dpLangevin->check_boundary_conditions(p))
    {
        std::cout << "Failure: wrong number of boundary conditions" << std::endl;
        return false;
    }
    if (not choose_integrator())
    { 
        std::cout << "Failure: unable to choose integrator" << std::endl;
        return false; 
    }        
    is_initialized = true;
    return is_initialized;
}

bool SimDP::run(const int n_next_epochs)
{
    if (not is_initialized) 
    { 
        std::cout << "Failure: must initialize first" << std::endl;
        return false; 
    }
    did_integrate = integrate(n_next_epochs);
    return did_integrate;
}

bool SimDP::postprocess()
{
    if (not is_initialized) 
    { 
        std::cout << "Failure: no data to process yet" << std::endl;
        return false; 
    }
    bool did_process = (
        prep_density_grid() and prep_t_epochs() and prep_mean_densities() 
    ); 
    return did_process;
}

int SimDP::get_n_epochs() const { return n_epochs; }
int SimDP::get_i_current_epoch() const { return i_current_epoch; }
int SimDP::get_i_next_epoch() const { return i_next_epoch; }
double SimDP::get_t_current_epoch() const { return t_current_epoch; }
double SimDP::get_t_next_epoch() const { return t_next_epoch; }
py_array_t SimDP::get_t_epochs() const { return pyarray_t_epochs; }
py_array_t SimDP::get_mean_densities() const { return pyarray_mean_densities; }
py_array_t SimDP::get_density() const { return pyarray_density; }

