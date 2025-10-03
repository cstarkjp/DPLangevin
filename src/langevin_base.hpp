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
class Langevin
{
protected:
    //! Runge-Kutta variable
    dbl_vec_t k1;
    //! Runge-Kutta variable
    dbl_vec_t k2;
    //! Runge-Kutta variable
    dbl_vec_t k3;
    //! Runge-Kutta variable
    dbl_vec_t k4;

    dbl_vec_t aux_cell_old, aux_cell_new;

    double dt, dx, dtm, dts;

    //! Dornic method stochastic-step variable
    double lambda;
    //! Dornic method stochastic-step variable
    double lambda_product;

    //! Function generating Poisson variates
    int_poisson_dist_t poisson;
    //! Function generating gamma variates
    dbl_gamma_dist_t gamma;
    //! Function generating normal variates
    dbl_normal_dist_t normal;

    //! Total number of cells in n-D grid
    int n_cells;
    //! Density field grid
    dbl_vec_t cell_density;
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
    Langevin() = default;
    //! Method to build 1d density grid and associated grid cell topology vectors
    bool construct_1D_grid(const Parameters parameters);
    //! Method to build 2d density grid and associated grid cell topology vectors
    bool construct_2D_grid(const Parameters parameters);
    //! Method to assign uniform random numbers to the initial density field
    void ic_random_uniform(
        rng_t &rng, const double min_value = 0.0, const double max_value = 1.0
    );
    void ic_constant_value(const double density_value=1.0);
    void ic_single_seed(const int i_node, const double value=1.0);
    void set_coefficients(const Coefficients &coefficients);
    void set_essential_coefficients(const Coefficients &coefficients);
    void set_lambdas();
    void integrate_rungekutta(rng_t &rng);
    void integrate_euler(rng_t &rng);
    void rk_f1(dbl_vec_t &aux_cell, dbl_vec_t &k1);
    void rk_f2f3(
        const dbl_vec_t &aux_old, 
        dbl_vec_t &aux_new, 
        dbl_vec_t &k_out, 
        const double dt_in
    );
    void rk_f4_and_stochastic(
        const dbl_vec_t &aux_old, 
        const dbl_vec_t &k1, 
        const dbl_vec_t &k2, 
        const dbl_vec_t &k3, 
        rng_t &rng
    );
    void euler_and_stochastic(dbl_vec_t &aux, rng_t &rng);
    double get_cell_density(const int) const;
    double get_mean_density() const;
    double get_poisson_mean() const;

    //! Method to set nonlinear coefficients for deterministic integration step: to be defined by application
    virtual void set_nonlinear_coefficients(const Coefficients &coefficients) {};
    //! Method to set nonlinear RHS of Langevin equation for deterministic integration step: to be defined by application
    virtual double nonlinear_rhs(const int i_cell, const dbl_vec_t &field) const 
        { return 0; };
};

#endif