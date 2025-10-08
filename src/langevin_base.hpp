/**
 * @file langevin_base.hpp
 * @brief Base class for Langevin equation integrator.
 */

#ifndef BASE_HPP
#define BASE_HPP

#include "general_coefficients.hpp"
#include "general_parameters.hpp"

/**
 * @brief Base class for Langevin equation integrator.
 * 
 * 
 */
class BaseLangevin
{
protected:
    //! Runge-Kutta variable #1
    dbl_vec_t k1;
    //! Runge-Kutta variable #2
    dbl_vec_t k2;
    //! Runge-Kutta variable #3
    dbl_vec_t k3;
    //! Runge-Kutta variable #4
    dbl_vec_t k4;

    //! Temporary density grid used to perform an integration step
    dbl_vec_t density_grid_aux_old;
    //! Temporary density grid used to perform an integration step
    dbl_vec_t density_grid_aux_new;

    //! Time step, i.e, epoch-to-epoch Δt
    double dt;
    //! Grid spacing, i.e., spacing Δx between cell centers in all directions
    double dx;
    //! Time step TBD
    double dtm;
    //! Time step TBD
    double dts;

    //! Dornic method stochastic-step variable
    double lambda;
    //! Dornic method stochastic-step variable
    double lambda_product;

    //! Function generating Poisson variates
    int_poisson_dist_t poisson_rng;
    //! Function generating gamma variates
    dbl_gamma_dist_t gamma_rng;
    //! Function generating normal variates
    dbl_normal_dist_t normal;

    //! Total number of cells in n-D grid
    int n_cells;
    //! Density field grid
    dbl_vec_t density_grid;
    //! Grid-average of density field
    double mean_density;

    //! Neighorhood topology for all grid cells
    std::vector< int_vec_t > neighbors;

    //! Dornic method coefficient
    double linear_coeff;
    //! Dornic method coefficient
    double noise_coeff;

public:
    //! Default constructor
    BaseLangevin() = default;
    //! Build 1d Langevin density field grid & topology
    bool construct_1D_grid(const Parameters parameters);
    //! Build 2d Langevin density field grid & topology
    bool construct_2D_grid(const Parameters parameters);
    //! Build 2d Langevin density field grid & mixedtopology
    bool construct_2D_grid_mixedtopology(const Parameters parameters);
    //! Initial condition for density field: uniformly random
    void ic_random_uniform(
        rng_t &rng, const double min_value = 0.0, const double max_value = 1.0
    );
    //! Initial condition for density field: uniformly constant
    void ic_constant_value(const double density_value=1.0);
    //! Initial condition for density field: single non-zero value
    void ic_single_seed(const int i_node, const double value=1.0);
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
    void rk_f1(dbl_vec_t &density_grid_aux, dbl_vec_t &k1);
    //! Parts #2 and #3 of Runge-Kutta integration step
    void rk_f2f3(
        const dbl_vec_t &density_grid_aux_old, 
        dbl_vec_t &density_grid_aux_new, 
        dbl_vec_t &k_out, 
        const double dt_in
    );
    //! Part #4 of Runge-Kutta integration step + stochastic step
    void rk_f4_and_stochastic(
        const dbl_vec_t &density_grid_aux_old, 
        const dbl_vec_t &k1, 
        const dbl_vec_t &k2, 
        const dbl_vec_t &k3, 
        rng_t &rng
    );
    //! Explicit Euler + stochastic integration
    void euler_and_stochastic(dbl_vec_t &density_grid_aux, rng_t &rng);
    //! Expose density grid at a given "index" (1d-flattened position of grid element)
    double get_density_grid_value(const int) const;
    //! Expose mean density
    double get_mean_density() const;
    //! Compute Poisson RNG mean
    double get_poisson_mean() const;

    //! Method to set nonlinear coefficients for deterministic integration step: to be defined by application
    virtual void set_nonlinear_coefficients(const Coefficients &coefficients) {};
    //! Method to set nonlinear RHS of Langevin equation for deterministic integration step: to be defined by application
    virtual double nonlinear_rhs(const int i_cell, const dbl_vec_t &field) const 
        { return 0; };
};

#endif