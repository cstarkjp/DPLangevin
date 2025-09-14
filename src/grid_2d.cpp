// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void Langevin::construct_2D_grid(const Parameters parameters)
{
    const int n_x = parameters.n_x;
    const int n_y = parameters.n_y;
    std::cout << "n_x: " << n_x << std::endl;
    std::cout << "n_y: " << n_y << std::endl;
    int i_cell, i_up, i_down, i_right, i_left;
    int i_top_row, i_bottom_row, i_right_column, i_left_column;

    neighbors = std::vector<int_vec_t>(n_x*n_y, int_vec_t(4));

    if (parameters.grid_topology==GridTopology::PERIODIC)
    {
        // Periodic grid topology in both x and y
        for (auto y=0; y<n_y; y++)
        {
            for (auto x=0; x<n_x; x++)
            {
                i_cell = x + y*n_x;
                i_up    = (y < n_y-1) ? i_cell + n_x : x;
                i_down  = (y > 0)     ? i_cell - n_x : x + (n_y-1)*n_x;
                i_right = (x < n_x-1) ? i_cell + 1 : 0 + y*n_x;
                i_left  = (x > 0)     ? i_cell - 1 : n_x-1 + y*n_x;

                neighbors[i_cell][0] = i_up;    // Up
                neighbors[i_cell][1] = i_down;  // Down
                neighbors[i_cell][2] = i_right; // Right   (VMB: left)
                neighbors[i_cell][3] = i_left;  // Left    (VMB: right)
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
            i_bottom_row = x;
            i_top_row    = x + (n_y-1)*n_x;

            // Each boundary cell has only 3 neighbors
            neighbors[i_bottom_row] = int_vec_t(3);
            neighbors[i_top_row]    = int_vec_t(3);

            neighbors[i_bottom_row][0] = i_bottom_row + n_x;
            neighbors[i_bottom_row][1] = i_bottom_row - 1; 
            neighbors[i_bottom_row][2] = i_bottom_row + 1; 
            neighbors[i_top_row][0] = i_top_row - n_x;
            neighbors[i_top_row][1] = i_top_row - 1; 
            neighbors[i_top_row][2] = i_top_row + 1;
        }
        // Left and right columns
        for (auto y=1; y<n_y-1; y++)
     {
            i_left_column  = y*n_x;
            i_right_column = (n_x-1) + y*n_x;

            // Each boundary cell has only 3 neighbors
            neighbors[i_left_column]  = int_vec_t(3);
            neighbors[i_right_column] = int_vec_t(3);

            neighbors[i_left_column][0] = i_left_column + n_x;
            neighbors[i_left_column][1] = i_left_column - n_x; 
            neighbors[i_left_column][2] = i_left_column + 1;
            neighbors[i_right_column][0] = i_right_column + n_x;
            neighbors[i_right_column][1] = i_right_column - n_x; 
            neighbors[i_right_column][2] = i_right_column - 1;
        }

        // Each corner cell has only 2 neighbors
        // Bottom-left corner
        i_cell = 0;
        neighbors[i_cell] = int_vec_t(2);
        neighbors[i_cell][0] = i_cell + 1; 
        neighbors[i_cell][1] = i_cell + n_x;
        // Bottom-right corner
        i_cell = n_x-1;
        neighbors[i_cell] = int_vec_t(2);
        neighbors[i_cell][0] = i_cell - 1;
        neighbors[i_cell][1] = i_cell + n_x;
        // Top-left corner
        i_cell = (n_y-1)*n_x;
        neighbors[i_cell] = int_vec_t(2);
        neighbors[i_cell][0] = i_cell + 1; 
        neighbors[i_cell][1] = i_cell - n_x;
        // Top-right corner
        i_cell = (n_x-1)+(n_y-1)*n_x;
        neighbors[i_cell] = int_vec_t(2);
        neighbors[i_cell][0] = i_cell - 1; 
        neighbors[i_cell][1] = i_cell - n_x;
    }
}