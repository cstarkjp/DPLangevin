/**
 * @file langevin_prepare.cpp
 * @brief Methods for setting the main Langevin model equation coefficients.
 */

#include "langevin_types.hpp"
#include "langevin_coefficients.hpp"
#include "langevin_base.hpp"

//! Set the Langevin equation coefficients and "lambda" coefficients
void BaseLangevin::prepare(const Coefficients& coefficients)
{
    linear_coefficient = coefficients.linear;
    noise_coefficient = coefficients.noise;
    const auto explcdt = exp(-linear_coefficient*dt);

    // Used in sampling gamma distribution
    lambda = (
        (2*linear_coefficient*explcdt)
            / 
        ((1-explcdt)*(noise_coefficient*noise_coefficient))
    );
    lambda_on_explcdt = lambda / explcdt;

    set_nonlinear_coefficients(coefficients);
}