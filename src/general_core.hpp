/**
 * @file general_core.hpp
 * @brief Core definitions for BaseLangevin integrator.
 */

#ifndef CORE_HPP
#define CORE_HPP

#include <iostream>
#include <vector>
#include <random>
#include <pybind11/numpy.h>
#include <valarray>

//! Use Mersenne Twister random number generator
typedef std::mt19937 rng_t;

//! Type for vectors of doubles
typedef std::vector<double> dbl_vec_t;
//! Type for vectors of integers
typedef std::vector<int> int_vec_t;

//! Type for density grid
typedef std::vector<double> grid_t;
// typedef std::valarray<double> grid_t;  // doesn't work

//! Type for density grid wiring
typedef std::vector<int> neighborhood_t;
//! Type for grid-cell neighborhood connections
typedef std::vector< neighborhood_t > grid_wiring_t;
// typedef std::valarray<int> grid_connection_t;  // doesn't work
// typedef std::valarray< grid_connection_t > grid_wiring_t;  // doesn't work

//! Type for function generating Poisson variates
typedef std::poisson_distribution<int> poisson_dist_t;
//! Type for function generating gamma variates
typedef std::gamma_distribution<double> gamma_dist_t;
//! Type for function generating Gaussian variates
typedef std::normal_distribution<double> normal_dist_t;
//! Type for function generating uniformly distributed variates
typedef std::uniform_real_distribution<double> uniform_dist_t;

#include "general_enums.hpp"
//! Type for specifying grid topology in each direction x, y, z...
typedef std::vector<GridTopology> gt_vec_t;
#include "general_coefficients.hpp"
#include "general_parameters.hpp"
#include "langevin_base.hpp"

namespace py = pybind11;
//! Type for Python arrays of doubles
typedef py::array_t<double, py::array::c_style> py_array_t;

//! Option to approximate Poisson distribution using Gaussian for large means
#ifndef APPROXIMATE_POISSON_DISTBN
#define APPROXIMATE_POISSON_DISTBN false
//! Threshold mean value above which Gaussian approximation should be used
#ifndef MU_THRESHOLD
#define MU_THRESHOLD 100.0
#endif
#endif

#endif