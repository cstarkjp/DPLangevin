// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef APPLICATIONS_HPP
#define APPLICATIONS_HPP

namespace py = pybind11;
typedef py::array_t<double, py::array::c_style> results_t;

auto dp_demo(
    double linear, double quadratic, double diffusion, double noise, 
    int n_cells, double t_max, double dx, double dt, int random_seed
) -> results_t;

#endif