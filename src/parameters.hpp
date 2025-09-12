// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef PARAMS_HPP
#define PARAMS_HPP
struct Parameters 
{
public:
    int n_cells;
    double t_max;
    double dx;
    double dt;
    int random_seed;
    GridDimension grid_dimension;
    GridTopology grid_topology;
    BoundaryCondition boundary_condition;
    InitialCondition initial_condition;
    IntegrationMethod integration_method;

    Parameters(
        int a, double b, double c, double d, int e, 
        GridDimension f, 
        GridTopology g,
        BoundaryCondition h,
        InitialCondition i,
        IntegrationMethod j
    ) : 
    n_cells(a), t_max(b), dx(c), dt(d), random_seed(e),
    grid_dimension(f), 
    grid_topology(g),
    boundary_condition(h),
    initial_condition(i), 
    integration_method(j)
    {}

    std::string gdstr(GridDimension gd) {
        switch (gd) {
            case GridDimension::D1: return "1d";
            case GridDimension::D2: return "2d";
            case GridDimension::D3: return "3d";
            default: return "Unknown";
        }
    }

    std::string gtstr(GridTopology gt) {
        switch (gt) {
            case GridTopology::BOUNDED: return "bounded";
            case GridTopology::PERIODIC: return "periodic";
            default: return "Unknown";
        }
    }

    std::string bcstr(BoundaryCondition bc) {
        switch (bc) {
            case BoundaryCondition::FLOATING: return "floating";
            case BoundaryCondition::FIXED_VALUE: return "fixed value";
            case BoundaryCondition::FIXED_FLUX: return "fixed flux";
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

    std::string icstr(IntegrationMethod ic) {
        switch (ic) {
            case IntegrationMethod::EULER: return "Euler";
            case IntegrationMethod::RUNGE_KUTTA: return "Runge-Kutta";
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
        std::cout<< "grid_topology: " 
            << gtstr(grid_topology) << std::endl;
        std::cout<< "boundary_condition: " 
            << bcstr(boundary_condition) << std::endl;
        std::cout<< "initial_condition: " 
            << icstr(initial_condition) << std::endl;
        std::cout<< "integration_method: " 
            << icstr(integration_method) << std::endl;
    }

    ~Parameters() {}
};

#endif