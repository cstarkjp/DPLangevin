/**
 * @file application_dplangevin.hpp
 * @brief DPLangevin model application of BaseLangevin class integrator.
 */

#ifndef DPLANGEVIN_HPP
#define DPLANGEVIN_HPP

/**
 * @brief DPLangevin model application of BaseLangevin class integrator.
 */
class DPLangevin : public BaseLangevin 
{
public:
    //! Coefficient in nonlinear term -bρ² in DP-Langevin equation
    double quadratic_coeff;
    //! Diffusion coefficient D in DP-Langevin equation
    double D;
    //! Constructor assuming default model parameters
    DPLangevin() = default;
    //! Constructor when model parameters are passed by the user
    DPLangevin(Parameters p);
    void set_nonlinear_coefficients(const Coefficients &coefficients);
    double nonlinear_rhs(const int i_cell, const dbl_vec_t &field) const;
};

#endif
