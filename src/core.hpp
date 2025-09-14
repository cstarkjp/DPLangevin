// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef CORE_HPP
#define CORE_HPP

#include <iostream>
#include <vector>
#include <random>
#include <pybind11/numpy.h>

// Option to approximate Poisson distribution using Gaussian for large means
#ifndef APPROXIMATE_POISSON_DISTBN
#define APPROXIMATE_POISSON_DISTBN false
#ifndef MU_THRESHOLD
#define MU_THRESHOLD 100.0
#endif
#endif

// Use Mersenne Twister random number generator
typedef std::mt19937 rng_t;
typedef std::vector<double> dbl_vec_t;
typedef std::vector<int> int_vec_t;
typedef std::poisson_distribution<int> int_poisson_dist_t;
typedef std::gamma_distribution<double> dbl_gamma_dist_t;
typedef std::normal_distribution<double> dbl_normal_dist_t;
typedef std::uniform_real_distribution<double> dbl_uniform_dist_t;

#include "enums.hpp"
#include "coefficients.hpp"
#include "parameters.hpp"
#include "langevin_base.hpp"

namespace py = pybind11;
typedef py::array_t<double, py::array::c_style> py_array_t;

#endif