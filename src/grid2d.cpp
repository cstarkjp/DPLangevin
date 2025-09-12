// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include "core.hpp"

void DornicBase::construct_2D_grid(const Parameters parameters)
{
    int top_row, bottom_row, right_column, left_column;
    const int n_x = parameters.grid_size.at(0);
    const int n_y = parameters.grid_size.at(1);
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
                auto i = x + y*n_x;
                top_row      = (y < n_y-1) ? x + (y+1)*n_x : x;
                bottom_row   = (y > 0)     ? x + (y-1)*n_x : x + (n_y-1)*n_x;
                right_column = (x < n_x-1) ? x+1 + y*n_x : 0 + y*n_x;
                left_column  = (x > 0)     ? x-1 + y*n_x : n_x-1 + y*n_x;

                neighbors[i][0] = top_row;    // Up
                neighbors[i][1] = bottom_row;  // Down
                neighbors[i][2] = right_column; // Left (CPS: huh??)
                neighbors[i][3] = left_column;  // Right (CPS: huh??)
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
                neighbors[i][0] = i + n_x;  // Up:   i + n_x // x + (y+1)*n_x;
                neighbors[i][1] = i - n_x;  // Down: i - n_x // x + (y-1)*n_x;
                neighbors[i][2] = i + 1;  // Right: i+1   (VMB: left)  // (x+1) + y*n_x;
                neighbors[i][3] = i - 1;  // Left:  i-1   (VMB: right) // (x-1) + y*n_x;
            }
        }
        // Top and bottom rows
        for (auto i=1; i<n_x-1; i++)
        {
            bottom_row = i;
            top_row    = i + (n_y-1)*n_x;

            // Each one of these have three neighbors, redefine
            neighbors[bottom_row] = int_vector(3);
            neighbors[top_row]    = int_vector(3);

            // Set the neighbors
            neighbors[bottom_row][0] = bottom_row + n_x;
            neighbors[bottom_row][1] = bottom_row - 1; 
            neighbors[bottom_row][2] = bottom_row + 1; 
            neighbors[top_row][0] = top_row - n_x;
            neighbors[top_row][1] = top_row - 1; 
            neighbors[top_row][2] = top_row + 1;
        }
        // Left and right columns
        for (auto i=1; i<n_y-1; i++)
        {
            left_column  = i*n_x;
            right_column = (n_x-1) + i*n_x;

            // Each one of these have three neighbors, redefine
            neighbors[left_column]  = int_vector(3);
            neighbors[right_column] = int_vector(3);

            // Set the neighbors
            neighbors[left_column][0] = left_column + n_x;
            neighbors[left_column][1] = left_column - n_x; 
            neighbors[left_column][2] = left_column + 1;
            neighbors[right_column][0] = right_column + n_x;
            neighbors[right_column][1] = right_column - n_x; 
            neighbors[right_column][2] = right_column - 1;
        }

        // Bottom-left corner
        neighbors[0][0] = 0 + 1; 
        neighbors[0][1] = 0 + n_x;
        // Bottom-right corner
        neighbors[n_x-1][0] = (n_x-1) - 1;
        neighbors[n_x-1][1] = (n_x-1) + n_x;
        // Top-left corner
        neighbors[(n_y-1)*n_x][0] = (n_y-1)*n_x + 1; 
        neighbors[(n_y-1)*n_x][1] = (n_y-2)*n_x;
        // Top-right corner
        neighbors[(n_x-1)+(n_y-1)*n_x][0] = (n_x-2)+(n_y-1)*n_x; 
        neighbors[(n_x-1)+(n_y-1)*n_x][1] = (n_x-1)+(n_y-2)*n_x;
    }
}