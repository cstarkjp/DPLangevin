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
    const double t_max, const double dx, const double dt, 
    const int random_seed,
    const GridDimension grid_dimension,
    const int_vec_t& grid_size,
    const GridTopology grid_topology,
    const BoundaryCondition boundary_condition,
    const InitialCondition initial_condition,
    const IntegrationMethod integration_method
) : f_coeffs(linear, quadratic, diffusion, noise),
    p( t_max, dx, dt, random_seed,
       grid_dimension, grid_size, grid_topology, 
       boundary_condition, initial_condition, integration_method )
{
    rng = new rng_t(p.random_seed); 
    dpLangevin = new DPLangevin(p);
    f_coeffs.print();
    p.print();
}

bool SimDP::initialize()
{
    if (not construct_grid()) { return false;}
    if (not initialize_grid()) { return false;}
    dpLangevin->set_coefficients(f_coeffs);
    is_initialized = true;
    return is_initialized;
}

bool SimDP::run(void)
{
    if (not is_initialized) 
    { 
        std::cout << "Failure: must initialize first" << std::endl;
        return false; 
    }
    n_epochs = count_epochs();
    epochs = dbl_vec_t(n_epochs, 0.0);
    mean_densities = dbl_vec_t(n_epochs, 0.0);
    bool did_integrate = integrate(epochs, mean_densities);
    bool did_finalize = (
        did_integrate and
        prep_epochs() and prep_mean_densities() and prep_density()
    ); 
    return did_finalize;
}

py_array_t SimDP::get_epochs() const { return return_epochs; }
py_array_t SimDP::get_mean_densities() const { return return_mean_densities; }
py_array_t SimDP::get_density() const { return return_density; }

