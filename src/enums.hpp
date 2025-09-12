// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef ENUMS_HPP
#define ENUMS_HPP

enum class GridDimension
{
    D1 = 1,
    D2 = 2,
    D3 = 3,
    D4 = 4
};

enum GridTopology
{
    BOUNDED = 1,
    PERIODIC = 2
};

enum BoundaryCondition
{
    FLOATING = 1,
    FIXED_VALUE = 2,
    FIXED_FLUX = 3
};

enum InitialCondition
{
    RANDOM_UNIFORM = 1,
    RANDOM_GAUSSIAN = 2,
    CONSTANT_VALUE = 3,
    SINGLE_SEED = 4
};

enum class IntegrationMethod
{
    EULER = 1,
    RUNGE_KUTTA = 2
};

#endif