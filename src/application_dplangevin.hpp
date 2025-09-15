// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef DPLANGEVIN_HPP
#define DPLANGEVIN_HPP

class DPLangevin : public Langevin 
{
public:
    double quadratic_coeff;
    double D;
    DPLangevin() = default;
    DPLangevin(Parameters p);
    void set_nonlinear_coefficients(const Coefficients &f_coefficients);
    double nonlinear_rhs(const int i_cell, const dbl_vec_t &field) const;
};

#endif
