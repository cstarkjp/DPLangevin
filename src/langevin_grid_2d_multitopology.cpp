/**
 * @file langevin_grid_2d.cpp
 * @brief Method for setting up a 2D grid for the model Langevin field.
 */

#include "general_core.hpp"

//! Construct 2D density field ρ(x,t) grid and corresponding cell-cell topologies
bool BaseLangevin::construct_2D_grid_multitopology(const Parameters p)
{
    const auto n_x = p.n_x;
    const auto n_y = p.n_y;
    int i_cell, i_top_row, i_bottom_row, i_right_column, i_left_column;

    // Flattened grid vector each with a set of 4 connection "i_node" elements.
    // Each i_node element will contain 1 of 4 possible neighbor locations.
    // Along edges these sets will be reduced to 3 elements.
    // At corners these sets will be reduced to 2 elements.
    neighbors = std::vector<int_vec_t>(n_x*n_y, int_vec_t(4));

    // Wiring lambdas
    auto connect_central_cell = [&](int x, int y)
    {
        // i_cell is the index of the flattened grid
        auto i_central_cell = x + y*n_x;
        // Each cell has 4 neighbors[i_cell] indexes
        neighbors[i_central_cell][0] = i_central_cell + n_x;  // Up:   i_cell + n_x // x + (y+1)*n_x;
        neighbors[i_central_cell][1] = i_central_cell - n_x;  // Down: i_cell - n_x // x + (y-1)*n_x;
        neighbors[i_central_cell][2] = i_central_cell + 1;    // Right: i_cell+1   (VMB: left)  // (x+1) + y*n_x;
        neighbors[i_central_cell][3] = i_central_cell - 1;    // Left:  i_cell-1   (VMB: right) // (x-1) + y*n_x;

    };
    auto connect_central_cells = [&]()
    {
        for (auto y=1; y<n_y-1; y++)
        {
            for (auto x=1; x<n_x-1; x++)
            {
                connect_central_cell(x,y);
            }
        }
    };
    auto connect_periodic_edge_cell_yplusminus = [&](int x, int y) 
    {
        auto i_edge_cell = x + y*n_x;
        auto i_yplus   = (y < n_y-1) ? i_edge_cell + n_x : x;
        auto i_yminus  = (y > 0)     ? i_edge_cell - n_x : x + (n_y-1)*n_x;
        neighbors[i_edge_cell][0] = i_yplus;   // Up
        neighbors[i_edge_cell][1] = i_yminus;  // Down
    };
    auto connect_periodic_edge_cell_xplusminus = [&](int x, int y) 
    {
        auto i_edge_cell = x + y*n_x;
        auto i_xplus  = (x < n_x-1) ? i_edge_cell + 1 : 0 + y*n_x;
        auto i_xminus = (x > 0)     ? i_edge_cell - 1 : n_x-1 + y*n_x;
        neighbors[i_edge_cell][2] = i_xplus;  // Right   (VMB: left)
        neighbors[i_edge_cell][3] = i_xminus; // Left    (VMB: right)
    };
    auto connect_periodic_edge_y_cells = [&](int x)
    {
        for (auto y=0; y<n_y; y++)
        {
            connect_periodic_edge_cell_yplusminus(x, y);
            connect_periodic_edge_cell_xplusminus(x, y);
        }
    };
    auto connect_periodic_edge_x_cells = [&](int y)
    {
        for (auto x=0; x<n_x; x++)
        {
            connect_periodic_edge_cell_yplusminus(x, y);
            connect_periodic_edge_cell_xplusminus(x, y);
        }
    };

    // Wiring of neighbors is simple everywhere except along the grid edges.
    // Here, do the wiring for all the non-edge grid cells.
    connect_central_cells();

    // Wire grid edge cells according to topology specs
    switch (p.grid_topology)
    {
        case GridTopology::PERIODIC:
            // Periodic grid topology in both x and y
            std::cout<< "construct_2D_grid_multitopology: " << "periodic" << std::endl;
            // Left edge, x=0, y cells
            connect_periodic_edge_y_cells(0);
            // Right edge, x=n_x-1, y cells
            connect_periodic_edge_y_cells(n_x-1);
            // Bottom edge, y=0, x cells
            connect_periodic_edge_x_cells(0);
            // Top edge, y=n_y-1, x cells
            connect_periodic_edge_x_cells(n_y-1);
            return true;

        case GridTopology::BOUNDED:
            // Bounded grid topology in both x and y
            std::cout<< "construct_2D_grid_multitopology: " << "bounded" << std::endl;
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
            return true;

        default:
            std::cout<< "construct_2D_grid_multitopology: " << "FAILED" << std::endl;
            return false;
    }
}