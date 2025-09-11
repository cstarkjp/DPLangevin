// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef STRUCT_HPP
#define STRUCT_HPP
// #include <string>

// Container for Langevin equation "coefficients"
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

enum class GridDimension
{
    D1 = 1,
    D2 = 2,
    D3 = 3,
    D4 = 4
};

enum InitialCondition
{
    RANDOM_UNIFORM = 1,
    RANDOM_GAUSSIAN = 2,
    CONSTANT_VALUE = 3,
    SINGLE_SEED = 4
};

enum BoundaryCondition
{
    PERIODIC = 1,
    FIXED_VALUE = 2,
    FIXED_FLUX = 3
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
    GridDimension grid_dimension;
    InitialCondition initial_condition;
    BoundaryCondition boundary_condition;

    Parameters(
        int a, double b, double c, double d, int e, 
        GridDimension f, 
        InitialCondition g,
        BoundaryCondition h
    ) : 
    n_cells(a), t_max(b), dx(c), dt(d), random_seed(e),
    grid_dimension(f), initial_condition(g), boundary_condition(h)
    {}

    std::string gdstr(GridDimension gd) {
        switch (gd) {
            case GridDimension::D1: return "1d";
            case GridDimension::D2: return "2d";
            case GridDimension::D3: return "3d";
            default: return "Unknown";
        }
    }

    std::string icstr(InitialCondition ic) {
        switch (ic) {
            case InitialCondition::RANDOM_UNIFORM: return "random uniform values";
            case InitialCondition::RANDOM_GAUSSIAN: return "random Gaussian values";
            case InitialCondition::CONSTANT_VALUE: return "constand value";
            case InitialCondition::SINGLE_SEED: return "single seed";
            default: return "Unknown";
        }
    }

    std::string bcstr(BoundaryCondition bc) {
        switch (bc) {
            case BoundaryCondition::PERIODIC: return "periodic";
            case BoundaryCondition::FIXED_VALUE: return "fixed value";
            case BoundaryCondition::FIXED_FLUX: return "fixed flux";
            default: return "Unknown";
        }
    }

    void print() {
        std::cout<< "n_cells: " << n_cells << std::endl;
        std::cout<< "t_max: " << t_max << std::endl;
        std::cout<< "dx: " << dx << std::endl;
        std::cout<< "dt: " << dt << std::endl;
        std::cout<< "random_seed: " << random_seed << std::endl;
        std::cout<< "grid_dimension: " 
            << gdstr(grid_dimension) << std::endl;
        std::cout<< "initial_condition: " 
            << icstr(initial_condition) << std::endl;
        std::cout<< "boundary_condition: " 
            << bcstr(boundary_condition) << std::endl;
    }

    ~Parameters() {}
};

#endif