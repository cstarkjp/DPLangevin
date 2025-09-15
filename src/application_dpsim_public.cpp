// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "general_core.hpp"
#include "application_dpsim.hpp"

SimDP::SimDP(
    const double linear, const double quadratic,
    const double diffusion, const double noise, 
    const double t_final, const double dx, const double dt, 
    const int random_seed,
    const GridDimension grid_dimension,
    const int_vec_t& grid_size,
    const GridTopology grid_topology,
    const BoundaryCondition boundary_condition,
    const InitialCondition initial_condition,
    const IntegrationMethod integration_method
) : coefficients(linear, quadratic, diffusion, noise),
    p( t_final, dx, dt, random_seed,
       grid_dimension, grid_size, grid_topology, 
       boundary_condition, initial_condition, integration_method )
{
    rng = new rng_t(p.random_seed); 
    dpLangevin = new DPLangevin(p);
    coefficients.print();
    p.print();
}

bool SimDP::initialize()
{
    if (not construct_grid()) { return false; }
    if (not initialize_grid()) { return false; }
    dpLangevin->set_coefficients(coefficients);
    n_epochs = count_epochs();
    epochs = dbl_vec_t(n_epochs, 0.0);
    mean_densities = dbl_vec_t(n_epochs, 0.0);
    // Treat the zeroth epoch as the initial grid state
    i_epoch = 1;
    t_epoch = p.dt;
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
    // Need to figure out how to segment sim right here
    // std::cout << "before i: " << i_epoch << std::endl;
    // std::cout << "before t: " << t_epoch << std::endl;
    did_integrate = integrate(n_next_epochs);
    // std::cout << "after  i: " << i_epoch << std::endl;
    // std::cout << "after  t: " << t_epoch << std::endl;
    return did_integrate;
}

bool SimDP::process()
{
    if (not is_initialized) 
    { 
        std::cout << "Failure: no data to process yet" << std::endl;
        return false; 
    }
    bool did_process = (
        prep_epochs() and prep_mean_densities() and prep_density()
    ); 
    return did_process;
}

int SimDP::get_n_epochs() const { return n_epochs; }
int SimDP::get_i_epoch() const { return i_epoch; }
double SimDP::get_t_epoch() const { return t_epoch; }
py_array_t SimDP::get_epochs() const { return return_epochs; }
py_array_t SimDP::get_mean_densities() const { return return_mean_densities; }
py_array_t SimDP::get_density() const { return return_density; }

