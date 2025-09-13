// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef LANGEVIN_BASE_HPP
#define LANGEVIN_BASE_HPP

#include "coefficients.hpp"
#include "parameters.hpp"

class LangevinBase
{
protected:
    // Runge-Kutta variables
    dbl_vec_t k1, k2, k3, k4, aux_cell_old, aux_cell_new;
    double dt, dx, dtm, dts;

    // Dornic method stochastic-step variables
    double lambda, lambda_product;
    int_poisson_dist_t poisson;
    dbl_gamma_dist_t gamma;
    dbl_normal_dist_t normal;

    // Grid variables
    int n_cells;
    dbl_vec_t cell_density;
    std::vector< int_vec_t > neighbors;
    double mean_density;

    // Dornic method coefficients
    double linear_coeff, noise_coeff;

public:
    LangevinBase(Parameters params)
    {
        dt = params.dt;
        dx = params.dx;
        dtm = 0.5*dt;
        dts = dt/6.0;
        n_cells = params.n_cells;
        cell_density = dbl_vec_t(n_cells, 0.0); 
        aux_cell_new = dbl_vec_t(n_cells);
        aux_cell_old = dbl_vec_t(n_cells);
        k1 = dbl_vec_t(n_cells, 0.0);
        k2 = dbl_vec_t(n_cells, 0.0);
        k3 = dbl_vec_t(n_cells, 0.0);
        k4 = dbl_vec_t(n_cells, 0.0);
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
        RNG &rng
    );
    void euler_and_stochastic(dbl_vec_t &aux, RNG &rng);
    double density(void);
    double avg_poisson_mean(void);

    // Defined by the application — these are placeholders
    virtual void set_nonlinear_coefficients(const Coefficients &f_coeffs) {};
    virtual auto nonlinear_rhs(const int i_cell, const dbl_vec_t &field) 
        const -> double{ return 0; };
};

#endif