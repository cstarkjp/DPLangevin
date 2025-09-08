// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

// Set all the parameters
void DornicBase::set_coefficients(const Coefficients &f_coeffs)
{
    set_essential_coefficients(f_coeffs);
    set_nonlinear_coefficients(f_coeffs);
    set_lambdas();
}

// Coefficients needed to integrate linear + noise parts
void DornicBase::set_essential_coefficients(const Coefficients &f_coeffs)
{
    linear_coeff = f_coeffs.linear;
    noise_coeff = f_coeffs.noise;
}

void DornicBase::set_lambdas(void)
{
    double lambda_const = 2.0/(noise_coeff*noise_coeff);
    double lambda_exp = exp(-linear_coeff*dt);
    lambda = lambda_const*linear_coeff*lambda_exp/(1.0-lambda_exp);
    lambda_product = lambda/lambda_exp;
}