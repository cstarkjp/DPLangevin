// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void DornicBase::construct_2D_grid(const Parameters parameters)
{
    int up, down, right, left;
    // const int L = sqrt(cell_density.size());
    const int n_x = parameters.grid_size.at(0);
    const int n_y = parameters.grid_size.at(1);
    const int L = n_x;
    std::cout << "n_x: " << n_x << std::endl;
    std::cout << "n_y: " << n_y << std::endl;

    neighbors = std::vector<int_vector>(n_cells, int_vector(4));

    if (parameters.grid_topology==GridTopology::PERIODIC)
    {
        // Periodic grid topology in both x and y
        for (auto y=0; y<n_y; y++)
        {
            for (auto x=0; x<n_x; x++)
            {
                up = (y < L-1) ? x + (y+1)*L : x;
                down = (y > 0) ? x + (y-1)*L : x + (L-1)*L;
                right = (x < L-1) ? x+1 + y*L : y*L;
                left = (x > 0) ? x-1 + y*L : L-1 + y*L;

                // up = (y < n_y-1) ? x + (y+1)*n_x : x;
                // down = (y > 0) ? x + (y-1)*n_x : x + (L-1)*n_x;
                // right = (x < L-1) ? x+1 + y*n_x : y*n_x;
                // left = (x > 0) ? x-1 + y*n_x : n_x-1 + y*n_x;

                auto i = x + y*n_x;
                neighbors[i][0] = up;    // Up
                neighbors[i][1] = down;  // Down
                neighbors[i][2] = right; // Left (CPS:??)
                neighbors[i][3] = left;  // Right (CPS:??)
            }
        }
    }
    else
    {
        // Bounded grid topology in both x and y
        for (auto y=1; y<n_y-1; y++)
        {
            for (auto x=1; x<n_x-1; x++)
            {
                // i is the index of the flattened grid
                auto i = x + y*n_x;
                // Each cell has 4 neighbors[i] indexes
                neighbors[i][0] = x + (y+1)*n_x;  // Up:   i + n_x
                neighbors[i][1] = x + (y-1)*n_x;  // Down: i - n_x
                neighbors[i][2] = (x+1) + y*n_x;  // Right: i+1   (VMB: left)
                neighbors[i][3] = (x-1) + y*n_x;  // Left:  i-1   (VMB: right)
            }
        }
        for (auto i=1; i<L-1; i++)
        {
            down  = i;             // Bottom row
            up    = i + (L-1)*L;   // Top row
            left  = i*L;           // Left column
            right = (L-1) + i*L;   // Right column

            // Each one of these have three neighbors, redefine
            neighbors[down]  = int_vector(3);
            neighbors[up]    = int_vector(3);
            neighbors[left]  = int_vector(3);
            neighbors[right] = int_vector(3);

            // Set the neighbors
            neighbors[down][0] = down + L;
            neighbors[down][1] = down - 1; 
            neighbors[down][2] = down + 1; 

            neighbors[up][0] = up - L;
            neighbors[up][1] = up - 1; 
            neighbors[up][2] = up + 1;

            neighbors[left][0] = left + L;
            neighbors[left][1] = left - L; 
            neighbors[left][2] = left + 1;

            neighbors[right][0] = right + L;
            neighbors[right][1] = right - L; 
            neighbors[right][2] = right - 1;
        }

        // Corners
        neighbors[0][0] = 1; 
        neighbors[0][1] = L;
        neighbors[L-1][0] = L-2; 
        neighbors[L-1][1] = (L-1)+L;
        neighbors[(L-1)*L][0] = 1+(L-1)*L; 
        neighbors[(L-1)*L][1] = (L-2)*L;
        neighbors[(L-1)+(L-1)*L][0] = (L-2)+(L-1)*L; 
        neighbors[(L-1)+(L-1)*L][1] = (L-1)+(L-2)*L;
    }
}