// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void LangevinBase::construct_2D_grid(const Parameters parameters)
{
    const int n_x = parameters.n_x;
    const int n_y = parameters.n_y;
    std::cout << "n_x: " << n_x << std::endl;
    std::cout << "n_y: " << n_y << std::endl;
    int i_cell, up_cell, down_cell, right_cell, left_cell;
    int top_row_cell, bottom_row_cell, right_column_cell, left_column_cell;

    neighbors = std::vector<int_vector>(n_x*n_y, int_vector(4));

    if (parameters.grid_topology==GridTopology::PERIODIC)
    {
        // Periodic grid topology in both x and y
        for (auto y=0; y<n_y; y++)
        {
            for (auto x=0; x<n_x; x++)
            {
                i_cell = x + y*n_x;
                up_cell    = (y < n_y-1) ? i_cell + n_x : x;
                down_cell  = (y > 0)     ? i_cell - n_x : x + (n_y-1)*n_x;
                right_cell = (x < n_x-1) ? i_cell + 1 : 0 + y*n_x;
                left_cell  = (x > 0)     ? i_cell - 1 : n_x-1 + y*n_x;

                neighbors[i_cell][0] = up_cell;    // Up
                neighbors[i_cell][1] = down_cell;  // Down
                neighbors[i_cell][2] = right_cell; // Right   (VMB: left)
                neighbors[i_cell][3] = left_cell;  // Left    (VMB: right)
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
                // i_cell is the index of the flattened grid
                i_cell = x + y*n_x;
                // Each cell has 4 neighbors[i_cell] indexes
                neighbors[i_cell][0] = i_cell + n_x;  // Up:   i_cell + n_x // x + (y+1)*n_x;
                neighbors[i_cell][1] = i_cell - n_x;  // Down: i_cell - n_x // x + (y-1)*n_x;
                neighbors[i_cell][2] = i_cell + 1;    // Right: i_cell+1   (VMB: left)  // (x+1) + y*n_x;
                neighbors[i_cell][3] = i_cell - 1;    // Left:  i_cell-1   (VMB: right) // (x-1) + y*n_x;
            }
        }
        // Top and bottom rows
        for (auto x=1; x<n_x-1; x++)
        {
            bottom_row_cell = x;
            top_row_cell    = x + (n_y-1)*n_x;

            // Each boundary cell has only 3 neighbors
            neighbors[bottom_row_cell] = int_vector(3);
            neighbors[top_row_cell]    = int_vector(3);

            neighbors[bottom_row_cell][0] = bottom_row_cell + n_x;
            neighbors[bottom_row_cell][1] = bottom_row_cell - 1; 
            neighbors[bottom_row_cell][2] = bottom_row_cell + 1; 
            neighbors[top_row_cell][0] = top_row_cell - n_x;
            neighbors[top_row_cell][1] = top_row_cell - 1; 
            neighbors[top_row_cell][2] = top_row_cell + 1;
        }
        // Left and right columns
        for (auto y=1; y<n_y-1; y++)
        {
            left_column_cell  = y*n_x;
            right_column_cell = (n_x-1) + y*n_x;

            // Each boundary cell has only 3 neighbors
            neighbors[left_column_cell]  = int_vector(3);
            neighbors[right_column_cell] = int_vector(3);

            neighbors[left_column_cell][0] = left_column_cell + n_x;
            neighbors[left_column_cell][1] = left_column_cell - n_x; 
            neighbors[left_column_cell][2] = left_column_cell + 1;
            neighbors[right_column_cell][0] = right_column_cell + n_x;
            neighbors[right_column_cell][1] = right_column_cell - n_x; 
            neighbors[right_column_cell][2] = right_column_cell - 1;
        }

        // Each corner cell has only 2 neighbors
        // Bottom-left corner
        i_cell = 0;
        neighbors[i_cell] = int_vector(2);
        neighbors[i_cell][0] = 0 + 1; 
        neighbors[i_cell][1] = 0 + n_x;
        // Bottom-right corner
        i_cell = n_x-1;
        neighbors[i_cell] = int_vector(2);
        neighbors[i_cell][0] = (n_x-1) - 1;
        neighbors[i_cell][1] = (n_x-1) + n_x;
        // Top-left corner
        i_cell = (n_y-1)*n_x;
        neighbors[i_cell] = int_vector(2);
        neighbors[i_cell][0] = (n_y-1)*n_x + 1; 
        neighbors[i_cell][1] = (n_y-2)*n_x;
        // Top-right corner
        i_cell = (n_x-1)+(n_y-1)*n_x;
        neighbors[i_cell] = int_vector(2);
        neighbors[i_cell][0] = (n_x-2)+(n_y-1)*n_x; 
        neighbors[i_cell][1] = (n_x-1)+(n_y-2)*n_x;
    }
}