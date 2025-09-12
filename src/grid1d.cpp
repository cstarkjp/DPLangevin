// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void DornicBase::construct_1D_grid(const Parameters parameters)
{
    const int n_x = parameters.grid_size.at(0);
    neighbors = std::vector<int_vector>(n_x, int_vector(2));

    // Everywhere except the grid ends
    for (auto i=1; i<n_x-1; i++)
    {
        // Each cell has a L and R neighbor whose indexes are specified here
        neighbors[i][0] = i-1;
        neighbors[i][1] = i+1;
    }

    // Grid ends
    if (parameters.grid_topology==GridTopology::PERIODIC)
    {
        // Each end cell neighbor is the other end cell, so wrap the indexes
        neighbors[0][0] = n_x-1;      // left-end left
        neighbors[0][1] = 1;          // left-end right
        neighbors[n_x-1][1] = n_x-2;  // right-end left VMB: [n_x-1][0] = n_x-2;
        neighbors[n_x-1][1] = 0;      // right-end right  
    }
    else
    {
        // Link each end cell to its adjacent cell only
        neighbors[0] = int_vector(1, 1);
        neighbors[n_x-1] = int_vector(n_x-2, n_x-2); // VMB: (1, n_x-2)
    }
}