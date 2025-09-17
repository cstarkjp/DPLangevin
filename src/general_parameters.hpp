// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP
struct Parameters 
{
public:
    const double t_final=0;
    const double dx=0;
    const double dt=0;
    const int random_seed=0;
    const GridDimension grid_dimension=GridDimension::D1;
    const int_vec_t& grid_size = {0};
    int n_cells=0;
    int n_x=0;
    int n_y=0;
    int n_z=0;  // but 3d, 4d grids not implemented yet (or ever)
    const GridTopology grid_topology=GridTopology::BOUNDED;
    const BoundaryCondition boundary_condition=BoundaryCondition::FLOATING;
    const InitialCondition initial_condition=InitialCondition::RANDOM_UNIFORM;
    const IntegrationMethod integration_method=IntegrationMethod::RUNGE_KUTTA;

    Parameters() = default;

    Parameters(
        const double t_final, 
        const double dx, const double dt, 
        const int rs, 
        const GridDimension gd, 
        const int_vec_t& gs,
        const GridTopology gt,
        const BoundaryCondition bc,
        const InitialCondition ic,
        const IntegrationMethod im
    ) : 
        t_final(t_final), 
        dx(dx), dt(dt), 
        random_seed(rs),
        grid_dimension(gd), 
        grid_size(gs),
        grid_topology(gt),
        boundary_condition(bc),
        initial_condition(ic), 
        integration_method(im)
    {
        n_x = gs.at(0);
        n_y = (gs.size()>1) ? gs.at(1) : 1;
        n_z = (gs.size()>2) ? gs.at(2) : 1;
        n_cells = n_x * n_y * n_z;
    }

    // Use overloading to provide alternate "printout" commands
    std::string report(GridDimension gd) 
    {
        switch (gd) {
            case GridDimension::D1: return "1d";
            case GridDimension::D2: return "2d";
            case GridDimension::D3: return "3d";
            default: return "Unknown";
        }
    }
    std::string report(GridTopology gt) 
    {
        switch (gt) {
            case GridTopology::BOUNDED: return "bounded";
            case GridTopology::PERIODIC: return "periodic";
            default: return "Unknown";
        }
    }
    std::string report(BoundaryCondition bc) 
    {
        switch (bc) {
            case BoundaryCondition::FLOATING: return "floating";
            case BoundaryCondition::FIXED_VALUE: return "fixed value";
            case BoundaryCondition::FIXED_FLUX: return "fixed flux";
            default: return "Unknown";
        }
    }
    std::string report(InitialCondition ic) 
    {
        switch (ic) {
            case InitialCondition::RANDOM_UNIFORM: return "random uniform values";
            case InitialCondition::RANDOM_GAUSSIAN: return "random Gaussian values";
            case InitialCondition::CONSTANT_VALUE: return "constant value";
            case InitialCondition::SINGLE_SEED: return "single seed";
            default: return "Unknown";
        }
    }
    std::string report(IntegrationMethod im) 
    {
        switch (im) {
            case IntegrationMethod::EULER: return "Euler";
            case IntegrationMethod::RUNGE_KUTTA: return "Runge-Kutta";
            default: return "Unknown";
        }
    }

    void print() 
    {
        std::cout<< "t_final: " << t_final << std::endl;
        std::cout<< "dx: " << dx << std::endl;
        std::cout<< "dt: " << dt << std::endl;
        std::cout<< "random_seed: " << random_seed << std::endl;
        std::cout<< "grid_dimension: " << report(grid_dimension) << std::endl;
        std::cout<< "grid_size: ";
        for (const auto& element : grid_size) {std::cout << element << " ";}
        std::cout<< std::endl;        
        std::cout<< "n_cells: " << n_cells << std::endl;
        std::cout<< "grid_topology: " << report(grid_topology) << std::endl;
        std::cout<< "boundary_condition: " << report(boundary_condition) << std::endl;
        std::cout<< "initial_condition: " << report(initial_condition) << std::endl;
        std::cout<< "integration_method: "  << report(integration_method) << std::endl;
    }

    ~Parameters() {}
};

#endif