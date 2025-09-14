// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef APPLICATION_DPSIM_HPP
#define APPLICATION_DPSIM_HPP

#include "application_dplangevin.hpp"

class SimDP 
{
private:
    Coefficients f_coeffs;
    Parameters p;
    RNG *rng; 
    DPLangevin *dpLangevin;
    int n_epochs;
    py_array_t return_epochs, return_mean_densities, density;
    dbl_vec_t epochs;
    dbl_vec_t mean_densities;
    bool is_initialized = false;
    void construct_grid();
    void initialize_grid();
    int count_epochs();
    bool integrate(dbl_vec_t& epochs, dbl_vec_t& mean_densities);
    bool prep_epochs();
    bool prep_mean_densities();
    bool prep_density();

public:
    SimDP(
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
    );
    bool initialize();
    bool run(void);
    py_array_t get_epochs() const;
    py_array_t get_mean_densities() const;
    py_array_t get_density() const;
};


#endif
