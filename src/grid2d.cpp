// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void DornicBase::construct_2D_grid(const Parameters parameters)
{
    // Only square grids are supported
    // TBD: implement rectangular grids
    const int L = sqrt(cell_density.size());
    int up, down, right, left;

    neighbors = std::vector<int_vector>(n_cells, int_vector(4));

    if (parameters.boundary_condition==BoundaryCondition::PERIODIC)
    {
        for (auto y=0; y<L; y++)
        {
            for (auto x=0; x<L; x++)
            {
                up = (y < L-1) ? x + (y+1)*L : x;
                down = (y > 0) ? x + (y-1)*L : x + (L-1)*L;
                right = (x < L-1) ? x+1 + y*L : y*L;
                left = (x > 0) ? x-1 + y*L : L-1 + y*L;

                auto i = x + y*L;
                neighbors[i][0] = up;    // Up
                neighbors[i][1] = down;  // Down
                neighbors[i][2] = right; // Left (CPS:??)
                neighbors[i][3] = left;  // Right (CPS:??)
            }
        }
    }
    else
    {
        for (auto y=1; y<L-1; y++)
        {
            for (auto x=1; x<L-1; x++)
            {
                auto i = x + y*L;
                neighbors[i][0] = x + (y+1)*L;  // Up
                neighbors[i][1] = x + (y-1)*L;  // Down
                neighbors[i][2] = (x+1) + y*L;  // Left (CPS:??)
                neighbors[i][3] = (x-1) + y*L;  // Right (CPS:??)
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