// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "general_core.hpp"

// Set all the parameters
void Langevin::set_coefficients(const Coefficients &coefficients)
{
    set_essential_coefficients(coefficients);
    set_nonlinear_coefficients(coefficients);
    set_lambdas();
}

// Coefficients needed to integrate linear + noise parts
void Langevin::set_essential_coefficients(const Coefficients &coefficients)
{
    linear_coeff = coefficients.linear;
    noise_coeff = coefficients.noise;
}

void Langevin::set_lambdas(void)
{
    double lambda_const = 2.0/(noise_coeff*noise_coeff);
    double lambda_exp = exp(-linear_coeff*dt);
    lambda = lambda_const*linear_coeff*lambda_exp/(1.0-lambda_exp);
    lambda_product = lambda/lambda_exp;
}