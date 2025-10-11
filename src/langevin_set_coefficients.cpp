/**
 * @file langevin_set_coefficients.cpp
 * @brief Methods for setting the main Langevin model equation coefficients.
 */

#include "general_types.hpp"

//! Set the Langevin equation coefficients and "lambda" coefficients
void BaseLangevin::set_coefficients(const Coefficients &coefficients)
{
    set_essential_coefficients(coefficients);
    set_nonlinear_coefficients(coefficients);
    set_lambdas();
}

//! Set the Langevin equation coefficients
void BaseLangevin::set_essential_coefficients(const Coefficients &coefficients)
{
    linear_coeff = coefficients.linear;
    noise_coeff = coefficients.noise;
}

//! Set "lambda" coefficients used in Poisson and gamma sampling
void BaseLangevin::set_lambdas(void)
{
    double lambda_const = 2.0/(noise_coeff*noise_coeff);
    double lambda_exp = exp(-linear_coeff*dt);
    // Used in sampling gamma distribution
    lambda = lambda_const*linear_coeff*lambda_exp/(1.0-lambda_exp);
    lambda_product = lambda/lambda_exp;
}