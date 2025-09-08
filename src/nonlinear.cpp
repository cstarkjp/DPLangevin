// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "base.hpp"

// Dummy definitions, to be replaced by application-specific versions

void Dornic::set_nonlinear_coefficients(const Coefficients &f_coefficients) {}
auto Dornic::nonlinear_rhs(const int i_node, const dbl_vector &field) 
    const -> double {}
