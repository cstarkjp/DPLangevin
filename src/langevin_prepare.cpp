/**
 * @file langevin_prepare.cpp
 * @brief Methods for setting the main Langevin model equation coefficients.
 */

#include "langevin_types.hpp"
#include "langevin_coefficients.hpp"
#include "langevin_base.hpp"

//! Set the Langevin equation coefficients and "lambda" coefficients
void BaseLangevin::prepare(const Coefficients &coefficients)
{
    linear_coeff = coefficients.linear;
    noise_coeff = coefficients.noise;

    double lambda_const = 2.0/(noise_coeff*noise_coeff);
    double lambda_exp = exp(-linear_coeff*dt);
    // Used in sampling gamma distribution
    lambda = lambda_const*linear_coeff*lambda_exp/(1.0-lambda_exp);
    lambda_product = lambda/lambda_exp;

    set_nonlinear_coefficients(coefficients);
}