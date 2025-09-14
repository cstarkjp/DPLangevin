// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef APPLICATION_DPLANGEVIN_HPP
#define APPLICATION_DPLANGEVIN_HPP

class DPLangevin : public LangevinBase 
{
public:
    double quadratic_coeff;
    double D;
    DPLangevin() = default;
    DPLangevin(Parameters params);
    void set_nonlinear_coefficients(const Coefficients &f_coefficients);
    auto nonlinear_rhs(const int i_cell, const dbl_vec_t &field) const ->double;
    void check(){std::cout << "DPLangevin::  check"  << std::endl;}
};

#endif
