/**
 * @file general_enums.hpp
 * @brief Enumerated parameter options.
 * 
 * Parameter options as enums, used in both C++ and Python, used to choose 
 * e.g. suitable grid geometries, topologies, boundary and initial conditions, 
 * and the lowest-level integration method.
 */

#ifndef ENUMS_HPP
#define ENUMS_HPP

enum class GridDimension
{
    D1 = 1,
    D2 = 2,
    D3 = 3,
    D4 = 4
};

enum class GridTopology
{
    BOUNDED = 1,
    PERIODIC = 2
};

enum class BoundaryCondition
{
    FLOATING = 1,
    FIXED_VALUE = 2,
    FIXED_FLUX = 3
};

enum class InitialCondition
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