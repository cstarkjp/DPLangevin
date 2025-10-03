/**
 * @file application_dplangevin.hpp
 * @brief Methods implementing deterministic, nonlinear part of Langevin eqn.
 */

#ifndef DPLANGEVIN_HPP
#define DPLANGEVIN_HPP

/**
 * @brief Methods implementing deterministic, nonlinear part of Langevin eqn.
 */
class DPLangevin : public Langevin 
{
public:
    double quadratic_coeff;
    double D;
    DPLangevin() = default;
    DPLangevin(Parameters p);
    void set_nonlinear_coefficients(const Coefficients &coefficients);
    double nonlinear_rhs(const int i_cell, const dbl_vec_t &field) const;
};

#endif
