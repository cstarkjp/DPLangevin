/**
 * @file application_dpsim_public.cpp
 * @brief Class to manage & execute DPLangevin model simulation: public methods.
 */ 

#include "general_core.hpp"
#include "application_dpsim.hpp"

SimDP::SimDP(
    const double linear, const double quadratic,
    const double diffusion, const double noise, 
    const double t_final, const double dx, const double dt, 
    const int random_seed,
    const GridDimension grid_dimension,
    const int_vec_t grid_size,
    const GridTopology grid_topology,
    const gt_vec_t grid_topologies,
    const BoundaryCondition boundary_condition,
    const InitialCondition initial_condition,
    const IntegrationMethod integration_method
) : coefficients(linear, quadratic, diffusion, noise),
    p(t_final, dx, dt, random_seed, 
        grid_dimension, grid_size, grid_topology, grid_topologies,
        boundary_condition, initial_condition, integration_method)
{
    rng = new rng_t(p.random_seed); 
    dpLangevin = new DPLangevin(p);
    coefficients.print();
    p.print();
    // std::cout 
    //     << "SimDP::SimDP: grid vec size: "
    //     << int(p.grid_size.size())
    //     << std::endl;
    // std::cout 
    //     << "SimDP::SimDP: gts vec size: "
    //     << int(p.grid_topologies.size())
    //     << std::endl;
}

bool SimDP::initialize()
{
    // std::cout 
    //     << "SimDP::initialize: grid: "
    //     << static_cast<int16_t>(this->p.grid_size.at(0))
    //     << ", "
    //     << static_cast<int16_t>(this->p.grid_size.at(1))
    //     << std::endl;
    // std::cout 
    //     << "SimDP::initialize: grid vec size: "
    //     << int(p.grid_size.size())
    //     << std::endl;
    // std::cout 
    //     << "SimDP::initialize: gts vec size: "
    //     << int(p.grid_topologies.size())
    //     << std::endl;
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
    // std::cout << "before i: " << i_next_epoch << std::endl;
    // std::cout << "before t: " << t_next_epoch << std::endl;
    did_integrate = integrate(n_next_epochs);
    // std::cout << "after  i: " << i_next_epoch << std::endl;
    // std::cout << "after  t: " << t_next_epoch << std::endl;
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
py_array_t SimDP::get_t_epochs() const { return return_t_epochs; }
py_array_t SimDP::get_mean_densities() const { return return_mean_densities; }
py_array_t SimDP::get_density() const { return return_density; }

