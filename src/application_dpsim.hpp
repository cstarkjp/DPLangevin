// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef SIMDP_HPP
#define SIMDP_HPP

#include "application_dplangevin.hpp"

class SimDP 
{
private:
    Coefficients coefficients;
    Parameters p;
    rng_t *rng; 
    DPLangevin *dpLangevin;
    int n_epochs;
    int i_epoch;
    double t_epoch;
    dbl_vec_t epochs, mean_densities;
    py_array_t return_t_epochs, return_mean_densities, return_density;
    bool did_integrate = false;
    bool is_initialized = false;
    bool construct_grid();
    bool initialize_grid();
    int count_epochs() const;
    bool integrate(const int n_next_epochs);
    bool prep_epochs();
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
    int get_i_epoch() const;
    double get_t_epoch() const;
    py_array_t get_t_epochs() const;
    py_array_t get_mean_densities() const;
    py_array_t get_density() const;
};


#endif
