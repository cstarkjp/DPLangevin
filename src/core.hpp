// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef BASE_HPP
#define BASE_HPP

#include <iostream>
#include <vector>
#include <random>

// Use Runge-Kutta integration method rather than Euler
#ifndef USE_RUNGEKUTTA
#define USE_RUNGEKUTTA true
#endif

// Option to approximate Poisson distribution using Gaussian for large means
#ifndef APPROXIMATE_POISSON_DISTBN
#define APPROXIMATE_POISSON_DISTBN false
#ifndef MU_THRESHOLD
#define MU_THRESHOLD 100.0
#endif
#endif

// Use Mersenne Twister random number generator
#ifndef RNG
#define RNG std::mt19937
#endif

typedef std::vector<double> dbl_vector;
typedef std::vector<int> int_vector;
typedef std::poisson_distribution<int> int_poisson_distbn;
typedef std::gamma_distribution<double> dbl_gamma_distbn;
typedef std::normal_distribution<double> dbl_normal_distbn;
typedef std::uniform_real_distribution<double> dbl_uniform_distbn;

#include "enums.hpp"
#include "coefficients.hpp"
#include "parameters.hpp"
#include "base.hpp"

#endif