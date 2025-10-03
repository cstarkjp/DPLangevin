/**
 * @file general_coefficients.hpp
 * @brief Container for nonlinear Langevin equation coefficients.
 */

#ifndef COEFFICIENTS_HPP
#define COEFFICIENTS_HPP

/**
 * @brief Container for nonlinear Langevin equation coefficients.
 *
 * Container for the set of coefficients in the nonlinear Langevin 
 * equation to be integrated, which here is the directed percolation (DP)
 * Langevin equation.
 * Includes a method to print out all coefficient values.
 *
 * @param linear Linear constant "a".
 * @param quadratic Nonlinear constant "b".
 * @param diffusion Diffusion rate.
 * @param noise Noise amplitude.
 * 
 */
struct Coefficients 
{
public:
    double linear;
    double quadratic;
    double diffusion;
    double noise;
    
    Coefficients(
        double a, double b, double c, double d
    ) : linear(a), quadratic(b), diffusion(c), noise(d) {}

    void print() {
        std::cout<< "linear: " << linear << std::endl;
        std::cout<< "quadratic: " << quadratic << std::endl;
        std::cout<< "diffusion: " << diffusion << std::endl;
        std::cout<< "noise: " << noise << std::endl;
    }

    ~Coefficients() {}
};

#endif