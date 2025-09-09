// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef STRUCT_HPP
#define STRUCT_HPP

// Container for Langevin equation "coefficients"
struct Coefficients 
{
public:
    double linear;
    double quadratic;
    double diffusion;
    double noise;

    Coefficients(double a, double b, double c, double d) 
        : linear(a), quadratic(b), diffusion(c), noise(d) {}

    void print() {
        std::cout<< "linear: " << linear << std::endl;
        std::cout<< "quadratic: " << quadratic << std::endl;
        std::cout<< "diffusion: " << diffusion << std::endl;
        std::cout<< "noise: " << noise << std::endl;
    }

    ~Coefficients() {}
};

// Container for model simulation parameters
struct Parameters 
{
public:
    int n_cells;
    double t_max;
    double dx;
    double dt;
    int random_seed;

    Parameters(int a, double b, double c, double d, int e) 
        : n_cells(a), t_max(b), dx(c), dt(d), random_seed(e) {}

    void print() {
        std::cout<< "n_cells: " << n_cells << std::endl;
        std::cout<< "t_max: " << t_max << std::endl;
        std::cout<< "dx: " << dx << std::endl;
        std::cout<< "dt: " << dt << std::endl;
        std::cout<< "random_seed: " << random_seed << std::endl;
    }

    ~Parameters() {}
};

#endif