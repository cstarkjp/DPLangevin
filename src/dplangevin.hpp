/**
 * @file dplangevin.hpp
 * @brief DPLangevin model application of BaseLangevin class integrator.
 */

#ifndef DPLANGEVIN_HPP
#define DPLANGEVIN_HPP

#include "langevin_types.hpp"
#include "langevin_base.hpp"

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
    //! Method to set nonlinear coefficients for deterministic integration step
    void set_nonlinear_coefficients(const Coefficients&  coefficients);
    //! Method to set nonlinear RHS of Langevin equation for deterministic integration step
    double nonlinear_rhs(const int i_cell, const grid_t&  field) const;
};

#endif
