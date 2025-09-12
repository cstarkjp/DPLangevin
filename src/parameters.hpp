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
    const double t_max;
    const double dx;
    const double dt;
    const int random_seed;
    const GridDimension grid_dimension;
    const int_vector& grid_size;
    int n_cells;
    int n_x;
    int n_y;
    int n_z;  // but 3d, 4d grids not implemented yet (or ever)
    const GridTopology grid_topology;
    const BoundaryCondition boundary_condition;
    const InitialCondition initial_condition;
    const IntegrationMethod integration_method;

    Parameters(
        const double b, const double c, const double d, const int e, 
        const GridDimension f, 
        const int_vector& k,
        const GridTopology g,
        const BoundaryCondition h,
        const InitialCondition i,
        const IntegrationMethod j
    ) : 
        t_max(b), dx(c), dt(d), random_seed(e),
        grid_dimension(f), 
        grid_size(k),
        grid_topology(g),
        boundary_condition(h),
        initial_condition(i), 
        integration_method(j)
    {
        n_x = k.at(0);
        n_y = (k.size()>1) ? k.at(1) : 1;
        n_z = (k.size()>2) ? k.at(2) : 1;
        n_cells = n_x * n_y * n_z;
    }

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
        std::cout<< "t_max: " << t_max << std::endl;
        std::cout<< "dx: " << dx << std::endl;
        std::cout<< "dt: " << dt << std::endl;
        std::cout<< "random_seed: " << random_seed << std::endl;
        std::cout<< "grid_dimension: " 
            << gdstr(grid_dimension) << std::endl;
        std::cout<< "n_cells: " << n_cells << std::endl;
        std::cout<< "grid_size: ";
        for (const auto& element : grid_size) {
            std::cout << element << " ";
        }
        std::cout<< std::endl;        
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