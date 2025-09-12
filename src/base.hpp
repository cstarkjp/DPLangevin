// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef DORNIC_HPP
#define DORNIC_HPP

#include "coefficients.hpp"
#include "parameters.hpp"

class DornicBase
{
protected:
    // Runge-Kutta variables
    dbl_vector k1, k2, k3, k4, aux_cell_old, aux_cell_new;
    double dt, dx, dtm, dts;

    // Dornic stochastic-step variables
    double lambda, lambda_product;
    int_poisson_distbn poisson;
    dbl_gamma_distbn gamma;
    dbl_normal_distbn normal;

    // Grid variables
    int n_cells;
    dbl_vector cell_density;
    std::vector< int_vector > neighbors;
    double mean_density;

    // Dornic method coefficients
    double linear_coeff, noise_coeff;

public:
    DornicBase(Parameters params)
    {
        // Initialization of integration increments
        dt = params.dt;
        dx = params.dx;

        // Initialization of Runge-Kutta variables
        dtm = 0.5*dt;
        dts = dt/6.0;

        // Initial number of cells, zero activity
        n_cells = params.n_cells;
        cell_density = dbl_vector(n_cells, 0.0); 
        aux_cell_new = dbl_vector(n_cells);
        aux_cell_old = dbl_vector(n_cells);

        // Assume constant cell numbers
        k1 = dbl_vector(n_cells, 0.0);
        k2 = dbl_vector(n_cells, 0.0);
        k3 = dbl_vector(n_cells, 0.0);
        k4 = dbl_vector(n_cells, 0.0);
    }

    void construct_1D_grid(const Parameters parameters);
    void construct_2D_grid(const Parameters parameters);
    void ic_random_uniform(
        RNG &rng, const double min_value = 0.0, const double max_value = 1.0
    );
    void ic_constant_value(const double density_value=1.0);
    void ic_single_seed(const int i_node, const double value=1.0);
    void set_coefficients(const Coefficients &f_coeffs);
    void set_essential_coefficients(const Coefficients &f_coeffs);
    void set_lambdas(void);
    void integrate_rungekutta(RNG &rng);
    void integrate_euler(RNG &rng);
    void rk_f1(dbl_vector &aux_cell, dbl_vector &k1);
    void rk_f2f3(
        const dbl_vector &aux_old, 
        dbl_vector &aux_new, 
        dbl_vector &k_out, 
        const double dt_in
    );
    void rk_f4_and_stochastic(
        const dbl_vector &aux_old, 
        const dbl_vector &k1, 
        const dbl_vector &k2, 
        const dbl_vector &k3, 
        RNG &rng
    );
    void euler_and_stochastic(dbl_vector &aux, RNG &rng);
    double density();
    double avg_poisson_mean();

    // Defined by the application — these are placeholders
    virtual void set_nonlinear_coefficients(const Coefficients &f_coeffs) {};
    virtual auto nonlinear_rhs(const int i_node, const dbl_vector &field) 
        const -> double{ return 0; };
};

#endif