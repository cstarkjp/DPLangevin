/**
 * @file langevin_grid_1d.cpp
 * @brief Method for setting up a 1D grid for the model Langevin field.
 */

#include "general_core.hpp"

bool Langevin::construct_1D_grid(const Parameters parameters)
{
    const int n_x = parameters.n_x;
    neighbors = std::vector<int_vec_t>(n_x, int_vec_t(2));

    // Everywhere except the grid ends
    for (auto i=1; i<n_x-1; i++)
    {
        // Each cell has a L and R neighbor whose indexes are specified here
        neighbors[i][0] = i-1;
        neighbors[i][1] = i+1;
    }

    // Grid ends
    switch (parameters.grid_topology)
    {
        case GridTopology::PERIODIC:
            // Each end cell neighbor is the other end cell, so wrap the indexes
            neighbors[0][0] = n_x-1;      // left-end left
            neighbors[0][1] = 1;          // left-end right
            neighbors[n_x-1][0] = n_x-2;  // right-end left VMB: [n_x-1][0] = n_x-2;
            neighbors[n_x-1][1] = 0;      // right-end right  
            
        case GridTopology::BOUNDED:
            // Link each end cell to its adjacent cell only
            neighbors[0] = int_vec_t(1, 1);
            neighbors[n_x-1] = int_vec_t(1, n_x-2);

        default:
            return false;
    }
}