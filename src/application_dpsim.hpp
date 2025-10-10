/**
 * @file application_dpsim.hpp
 * @brief Class to manage & execute model simulation using the DPLangevin integrator.
 */ 

#ifndef SIMDP_HPP
#define SIMDP_HPP

#include "application_dplangevin.hpp"

/**
 * @brief Class to manage & execute model simulation using the DPLangevin integrator.
 *
 * Class to manage & execute model simulation using an instances of 
 * the DPLangevin integrator class, the Coefficients struct, and the 
 * Parameters struct.
 */
class SimDP 
{
private:
    //! Langevin equation coefficients
    Coefficients coefficients;
    //! Model simulation parameters
    Parameters p;
    //! Random number generation function (Mersenne prime) (pointer to RNG)
    rng_t *rng; 
    //! Instance of DP Langevin integrator class (pointer to instance)
    DPLangevin *dpLangevin;
    //! Integrator: either a Runge-Kutta or an Euler method
    void (DPLangevin::*integrator)(rng_t&);
    //! Total number of simulation time steps aka "epochs"
    int n_epochs;
    //! Index of current epoch aka time step
    int i_current_epoch;
    //! Index of next epoch aka time step
    int i_next_epoch;
    //! Time of current epoch
    double t_current_epoch;
    //! Time of next epoch
    double t_next_epoch;
    //! Vector time-series of epochs
    dbl_vec_t t_epochs;
    //! Vector time-series of grid-averaged field density values
    dbl_vec_t mean_densities;
    //! Python-compatible array of epochs time-series
    py_array_t pyarray_t_epochs;
    //! Python-compatible array of mean density time-series
    py_array_t pyarray_mean_densities;
    //! Python-compatible array of current density grid
    py_array_t pyarray_density;
    //! Flag whether integration step was successful or not
    bool did_integrate = false;
    //! Flag whether simulation has been initialized or not
    bool is_initialized = false;

    //! Construct Langevin density field grid of appropriate n-D dimension
    bool construct_grid();
    //! Set initial condition of Langevin density field grid
    bool initialize_grid();
    //! Count upcoming number of epochs by running a dummy time-stepping loop
    int count_epochs() const;
    //! Choose integrator function implementing RK or Euler
    bool choose_integrator();
    //! Apply boundary conditions along edges
    void apply_boundary_conditions();
    //! Perform Dornic-type integration of the DP Langevin equation for `n_next_epochs`
    bool integrate(const int n_next_epochs);
    //! Generate a Python-compatible version of the epochs time-series vector
    bool prep_t_epochs();
    //! Generate a Python-compatible version of the mean densities time-series vector
    bool prep_mean_densities();
    //! Generate a Python-compatible version of the current density grid
    bool prep_density_grid();

public:
    //! Constructor
    SimDP(
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
    );
    //! Initialize the model simulation
    bool initialize();
    //! Execute the model simulation for `n_next_epochs`
    bool run(const int n_next_epochs);
    //! Process the model results data if available
    bool postprocess();
    //! Fetch the total number of simulation epochs
    int get_n_epochs() const;
    //! Fetch the index of the current epoch of the simulation
    int get_i_current_epoch() const;
    //! Fetch the index of the next epoch of the simulation
    int get_i_next_epoch() const;
    //! Fetch the current epoch (time) of the simulation
    double get_t_current_epoch() const;
    //! Fetch the next epoch (time) of the simulation
    double get_t_next_epoch() const;
    //! Fetch a times-series vector of the simulation epochs as a Python array
    py_array_t get_t_epochs() const;
    //! Fetch a times-series vector of the grid-averaged density field over time as a Python array
    py_array_t get_mean_densities() const;
    //! Fetch the current Langevin density field grid as a Python array
    py_array_t get_density() const;
};


#endif
