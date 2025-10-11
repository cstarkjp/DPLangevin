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
    grid_t density_grid_aux1;
    //! Temporary density grid used to perform an integration step
    grid_t density_grid_aux2;

public:
    //! Default constructor
    BaseLangevin() = default;
    //! Construct Langevin density field grid of appropriate n-D dimension
    bool construct_grid(const Parameters parameters);
    //! Build 1d Langevin density field grid & topology
    bool construct_1D_grid(const Parameters parameters);
    //! Build 2d Langevin density field grid & mixed topology
    bool construct_2D_grid(const Parameters parameters);
    //! Initial condition for density field: uniformly random
    void ic_random_uniform(
        rng_t& rng, 
        const double min_value=0.0, 
        const double max_value=1.0
    );
    //! Set initial condition of Langevin density field grid
    bool initialize_grid(const Parameters parameters, rng_t& rng);
    //! Initial condition for density field: uniformly constant
    void ic_constant_value(const double density_value=1.0);
    //! Initial condition for density field: single non-zero value
    void ic_single_seed(const int i, const double value=1.0);
    //! Method to set Langevin equation coefficients and "lambda" constants
    void set_coefficients(const Coefficients& coefficients);
    //! Method to set Langevin equation coefficients
    void set_essential_coefficients(const Coefficients& coefficients);
    //! Method to set "lambda" constants
    void set_lambdas();
    //! Check we have 2N boundary conditions for an N-dimensional grid
    bool check_boundary_conditions(const Parameters parameters);
    //! Set density field values only the grid edges per bc specs
    void apply_boundary_conditions(const Parameters parameters);
    //! Runge-Kutta + stochastic integration + grid update
    void integrate_rungekutta(rng_t& rng);
    //! Explicit Euler + stochastic integration + grid update
    void integrate_euler(rng_t& rng);
    //! Step #1 of Runge-Kutta integration
    void rk_f1(
        grid_t& density_grid_aux, 
        grid_t& k1_grid, 
        const double dtm
    );
    //! Steps #2 and #3 of Runge-Kutta integration
    void rk_f2f3(
        const grid_t& density_grid_aux_in, 
        grid_t& density_grid_aux_out, 
        grid_t& k_out, 
        const double dtm
    );
    //! Step #4 of Runge-Kutta integration + stochastic step
    void rk_f4_and_stochastic(
        const grid_t& density_grid_aux, 
        const grid_t& k1_grid, 
        const grid_t& k2_grid, 
        const grid_t& k3_grid, 
        rng_t& rng, 
        const double dtm
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