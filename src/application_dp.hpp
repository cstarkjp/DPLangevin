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

enum GridDimension
{
    D1 = 1,
    D2 = 2,
    D3 = 3,
    D4 = 4
};

enum Initialization
{
    RANDOM_UNIFORM,
    RANDOM_GAUSSIAN,
    CONSTANT_VALUE,
    SINGLE_SEED
};

enum BoundaryCondition
{
    PERIODIC,
    FIXED_VALUE,
    FIXED_FLUX
};

auto dp(
    double linear, double quadratic, double diffusion, double noise, 
    int n_cells, double t_max, double dx, double dt, int random_seed
) -> results_t;

#endif