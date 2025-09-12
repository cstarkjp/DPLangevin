// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void DornicBase::construct_1D_grid(const Parameters parameters)
{
    neighbors = std::vector<int_vector>(n_cells, int_vector(2));

    // Everywhere except the grid ends
    for (auto i=1; i<n_cells-1; i++)
    {
        neighbors[i][0] = i-1;
        neighbors[i][1] = i+1;
    }

    // Grid ends
    if (parameters.boundary_condition==BoundaryCondition::PERIODIC)
    {
        neighbors[0][0] = n_cells-1;
        neighbors[0][1] = 1;
        neighbors[n_cells-1][1] = n_cells-2;
        neighbors[n_cells-1][1] = 0;        
    }
    else
    {
        // These only have size one, redefine
        neighbors[0] = int_vector(1, 1);
        neighbors[n_cells-1] = int_vector(1, n_cells-2);
    }
}