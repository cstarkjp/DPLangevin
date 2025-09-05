// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef CLASS_HPP
#define CLASS_HPP

#include "struct.hpp"
class Dornic
{

private:
    // Runge-Kutta variables
    dbl_vector k1, k2, k3, k4;
    dbl_vector aux_cell_old, aux_cell_new;
    double dt, dx, dtm, dts;

    // Dornic stochastic step variables
    double lambda, lambda_product;
    int_poisson_distbn poisson;
    dbl_gamma_distbn gamma;
    dbl_normal_distbn normal;

    // Grid/network variables
    dbl_vector cell_density;
    double mean_density;
    int n_cells;
    std::vector< int_vector > neighbors;

    // Dornic method coefficients
    double linear_coeff, noise_coeff;

public:
    Dornic(Parameters params)
    {
        // Initialization of integration increments
        dt = params.dt; //dt_in;
        dx = params.dx; //dx_in;

        // Initialization of Runge-Kutta variables
        dtm = 0.5*dt;
        dts = dt/6.0;

        // Initial number of cells, set no activity for the system
        n_cells = params.n_cells; //total_n_cells;
        cell_density = dbl_vector(n_cells, 0.0); 
        aux_cell_new = dbl_vector(n_cells);
        aux_cell_old = dbl_vector(n_cells);

        // Assume constant cell numbers
        k1 = dbl_vector(n_cells, 0.0);
        k2 = dbl_vector(n_cells, 0.0);
        k3 = dbl_vector(n_cells, 0.0);
        k4 = dbl_vector(n_cells, 0.0);
    }

    // Declarations
    void integration(RNG &rng);
    void integration(RNG &rng, const Coefficients &f_coeffs);
    void integrate(RNG &rng);
    void set_coefficients(const Coefficients &f_coeffs);
    void set_essential_coefficients(const Coefficients &f_coeffs);
    void set_lambdas(void);
    void set_nonlinear_coefficients(const Coefficients &f_coeffs);
    void f1(dbl_vector &aux_cell, dbl_vector &k1);
    void f2f3(
        const dbl_vector &aux_old, 
        dbl_vector &aux_new, 
        dbl_vector &k_out, 
        const double dt_in
    );
    void f4_and_stochastic(
        const dbl_vector &aux_old, 
        const dbl_vector &k1, 
        const dbl_vector &k2, 
        const dbl_vector &k3, 
        RNG &rng
    );
    void euler_and_stochastic(dbl_vector &aux, RNG &rng);
    double nonlinear_rhs(const int i_node, const dbl_vector &field) const;
    void construct_1D_grid(const bool periodic = false);
    void construct_2D_grid(const bool periodic = false);
    void construct_custom_network(const std::vector<int_vector> &network);
    void random_intial_condition(
        RNG &rng, const double min_value = 0.0, const double max_value = 1.0
    );
    void homogeneous_initial_condition(const double density_value);
    void single_seed(const int i_node, const double value);
    double density();
    double avg_poisson_mean();
};

#endif