/**
 * @file langevin_base.hpp
 * @brief Base class for Langevin equation integrator.
 */

#ifndef BASE_HPP
#define BASE_HPP

#include "general_coefficients.hpp"
#include "general_parameters.hpp"

/**
 * @details Base class for Langevin equation integrator.
 */
class BaseLangevin
{
protected:
    //! Total number of cells in n-D grid
    int n_cells;
    //! Density field grid
    grid_t density_grid;
     //! Neighorhood topology for all grid cells
    grid_wiring_t grid_wiring;
   
    //! Time step, i.e, epoch-to-epoch Δt
    double dt;
    //! Grid spacing, i.e., spacing Δx between cell centers in all directions
    double dx;
    //! Grid-average of density field
    double mean_density;

    //! Function generating Poisson variates
    poisson_dist_t poisson_rng;
    //! Function generating gamma variates
    gamma_dist_t gamma_rng;
    //! Function generating normal variates
    normal_dist_t normal;

    //! Dornic method coefficient
    double linear_coeff;
    //! Dornic method coefficient
    double noise_coeff;
    //! Dornic method stochastic-step variable
    double lambda;
    //! Dornic method stochastic-step variable
    double lambda_product;

    //! Runge-Kutta variable grid #1
    grid_t k1_grid;
    //! Runge-Kutta variable grid #2
    grid_t k2_grid;
    //! Runge-Kutta variable grid #3
    grid_t k3_grid;
    //! Temporary density grid used to perform an integration step
    grid_t density_grid_aux_old;
    //! Temporary density grid used to perform an integration step
    grid_t density_grid_aux_new;

public:
    //! Default constructor
    BaseLangevin() = default;
    //! Build 1d Langevin density field grid & topology
    bool construct_1D_grid(const Parameters parameters);
    //! Build 2d Langevin density field grid & mixed topology
    bool construct_2D_grid(const Parameters parameters);
    //! Initial condition for density field: uniformly random
    void ic_random_uniform(
        rng_t &rng, 
        const double min_value=0.0, 
        const double max_value=1.0
    );
    //! Initial condition for density field: uniformly constant
    void ic_constant_value(const double density_value=1.0);
    //! Initial condition for density field: single non-zero value
    void ic_single_seed(const int i, const double value=1.0);
    //! Method to set Langevin equation coefficients and "lambda" constants
    void set_coefficients(const Coefficients &coefficients);
    //! Method to set Langevin equation coefficients
    void set_essential_coefficients(const Coefficients &coefficients);
    //! Method to set "lambda" constants
    void set_lambdas();
    //! Runge-Kutta + stochastic integration + grid update
    void integrate_rungekutta(rng_t &rng);
    //! Explicit Euler + stochastic integration + grid update
    void integrate_euler(rng_t &rng);
    //! Part #1 of Runge-Kutta integration step
    void rk_f1(grid_t &density_grid_aux, grid_t &k1_grid);
    //! Parts #2 and #3 of Runge-Kutta integration step
    void rk_f2f3(
        const grid_t &density_grid_aux_old, 
        grid_t &density_grid_aux_new, 
        grid_t &k_out, 
        const double dt_in
    );
    //! Part #4 of Runge-Kutta integration step + stochastic step
    void rk_f4_and_stochastic(
        const grid_t &density_grid_aux_old, 
        const grid_t &k1_grid, 
        const grid_t &k2_grid, 
        const grid_t &k3_grid, 
        rng_t &rng
    );
    //! Explicit Euler + stochastic integration
    void euler_and_stochastic(grid_t &density_grid_aux, rng_t &rng);
    //! Expose density grid at a given "index" (1d-flattened position of grid element)
    double get_density_grid_value(const int) const;
    //! Expose mean density
    double get_mean_density() const;
    //! Compute Poisson RNG mean
    double get_poisson_mean() const;

    //! Method to set nonlinear coefficients for deterministic integration step: to be defined by application
    virtual void set_nonlinear_coefficients(const Coefficients &coefficients) {};
    //! Method to set nonlinear RHS of Langevin equation for deterministic integration step: to be defined by application
    virtual double nonlinear_rhs(const int i_cell, const grid_t &field) const 
        { return 0; };
};

#endif