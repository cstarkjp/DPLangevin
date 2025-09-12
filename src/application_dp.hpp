// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef APPLICATION_DP_HPP
#define APPLICATION_DP_HPP

namespace py = pybind11;
typedef py::array_t<double, py::array::c_style> results_t;

results_t dp(
    double linear, double quadratic, double diffusion, double noise, 
    int n_cells, double t_max, double dx, double dt, int random_seed,
    GridDimension grid_dimension, 
    InitialCondition initial_condition,
    BoundaryCondition boundary_condition,
    IntegrationMethod integration_method
);

#endif