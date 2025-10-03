/**
 * @file application_dpsim.hpp
 * @brief Class to manage & execute Langevin model simulation.
 */ 

#ifndef SIMDP_HPP
#define SIMDP_HPP

#include "application_dplangevin.hpp"

/**
 * @brief Class to manage & execute Langevin model simulation.
 *
 * Class to manage & execute Langevin model simulation using DPLangevin integrator.
 */
class SimDP 
{
private:
    Coefficients coefficients;
    Parameters p;
    rng_t *rng; 
    DPLangevin *dpLangevin;
    int n_epochs;
    int i_next_epoch, i_current_epoch;
    double t_next_epoch, t_current_epoch;
    dbl_vec_t t_epochs, mean_densities;
    py_array_t return_t_epochs, return_mean_densities, return_density;
    bool did_integrate = false;
    bool is_initialized = false;
    bool construct_grid();
    bool initialize_grid();
    int count_epochs() const;
    bool integrate(const int n_next_epochs);
    bool prep_t_epochs();
    bool prep_mean_densities();
    bool prep_density();

public:
    SimDP(
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
    );
    bool initialize();
    bool run(const int n_next_epochs);
    bool process();
    int get_n_epochs() const;
    int get_i_current_epoch() const;
    int get_i_next_epoch() const;
    double get_t_current_epoch() const;
    double get_t_next_epoch() const;
    py_array_t get_t_epochs() const;
    py_array_t get_mean_densities() const;
    py_array_t get_density() const;
};


#endif
