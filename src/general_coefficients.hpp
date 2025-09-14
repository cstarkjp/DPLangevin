// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef COEFFS_HPP
#define COEFFS_HPP

struct Coefficients 
{
public:
    double linear;
    double quadratic;
    double diffusion;
    double noise;
    
    Coefficients(
        double a, double b, double c, double d
    ) : 
    linear(a), quadratic(b), diffusion(c), noise(d)
    {}

    void print() {
        std::cout<< "linear: " << linear << std::endl;
        std::cout<< "quadratic: " << quadratic << std::endl;
        std::cout<< "diffusion: " << diffusion << std::endl;
        std::cout<< "noise: " << noise << std::endl;
    }

    ~Coefficients() {}
};

#endif