/**
 * @file langevin_grid_2d.cpp
 * @brief Method for setting up a 2D grid for the model Langevin field.
 */

#include "general_core.hpp"

//! Construct 2D density field ρ(x,t) grid and corresponding cell-cell topologies
bool BaseLangevin::construct_2D_grid_multitopology(const Parameters p)
{
    // Shorthand
    const auto n_x = p.n_x;
    const auto n_y = p.n_y;

    // Flattened grid vector each with a set of <=4 connection i_node elements.
    // Each i_node element will link to 1 of <=4 possible neighbor locations.
    // Along edges these sets will be reduced to 3 elements.
    // At corners these sets will be reduced to 2 elements.
    neighbors = std::vector<int_vec_t>(n_x*n_y, int_vec_t(4));

    // Central cells
    // Single cell
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
    // All grid cells except edges
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
    
    // Periodic
    auto connect_periodic_edge_cell_yplus = [&](int x, int y) 
    {
        auto i_edge_cell = x + y*n_x;
        auto i_yplus = (y < n_y-1) ? i_edge_cell + n_x : x;
        neighbors[i_edge_cell][0] = i_yplus;   // Up
    };
    auto connect_periodic_edge_cell_yminus = [&](int x, int y) 
    {
        auto i_edge_cell = x + y*n_x;
        auto i_yminus = (y > 0) ? i_edge_cell - n_x : x + (n_y-1)*n_x;
        neighbors[i_edge_cell][1] = i_yminus;  // Down
    };
    auto connect_periodic_edge_cell_xplus = [&](int x, int y) 
    {
        auto i_edge_cell = x + y*n_x;
        auto i_xplus = (x < n_x-1) ? i_edge_cell + 1 : 0 + y*n_x;
        neighbors[i_edge_cell][2] = i_xplus;  // Right   (VMB: left)
    };
    auto connect_periodic_edge_cell_xminus = [&](int x, int y) 
    {
        auto i_edge_cell = x + y*n_x;
        auto i_xminus = (x > 0) ? i_edge_cell - 1 : n_x-1 + y*n_x;
        neighbors[i_edge_cell][3] = i_xminus; // Left    (VMB: right)
    };
    auto connect_periodic_edge_cells = [&](int x, int y) 
    {
        connect_periodic_edge_cell_yplus(x, y);
        connect_periodic_edge_cell_yminus(x, y);
        connect_periodic_edge_cell_xplus(x, y);
        connect_periodic_edge_cell_xminus(x, y);        
    };
    // x=left edge or x=right edge, loop over y cells
    auto connect_periodic_y_edge_cells = [&](int x)
    {
        assert(x==0 or x==n_x-1);
        for (auto y=1; y<n_y-1; y++)
        {
            connect_periodic_edge_cells(x, y);
        }
    };
    // y=bottom edge or y=top edge, loop over x cells
    auto connect_periodic_x_edge_cells = [&](int y)
    {
        assert(y==0 or y==n_y-1);
        for (auto x=1; x<n_x-1; x++)
        {
            connect_periodic_edge_cells(x, y);
        }
    };
    // Bottom-left corner
    auto connect_periodic_bl_corner = [&]()
    {
        connect_periodic_edge_cells(0, 0);
    };
    // Bottom-right corner
    auto connect_periodic_br_corner = [&]()
    {
        connect_periodic_edge_cells(n_x-1, 0);
    };
    // Top-left corner
    auto connect_periodic_tl_corner = [&]()
    {
        connect_periodic_edge_cells(0, n_y-1);
    };
    // Top-right corner
    auto connect_periodic_tr_corner = [&]()
    {
        connect_periodic_edge_cells(n_x-1, n_y-1);
    };


    // Bounded
    // Left and right edges, loop over y cells
    auto connect_bounded_y_edge_cells = [&](int x)
    {
        assert(x==0 or x==n_x-1);
        auto plus_or_minus = (x==0) ? +1 : -1;
        for (auto y=1; y<n_y-1; y++)
        {
            auto i_cell = x + y*n_x;
            neighbors[i_cell] = int_vec_t(3);
            neighbors[i_cell][0] = i_cell + n_x*plus_or_minus;
            neighbors[i_cell][1] = i_cell - n_x*plus_or_minus; 
            neighbors[i_cell][2] = i_cell + plus_or_minus;
        }
    };
    // Bottom and top edges, loop over x cells
    auto connect_bounded_x_edge_cells = [&](int y)
    {
        assert(y==0 or y==n_y-1);
        auto plus_or_minus = (y==0) ? +1 : -1;
        for (auto x=1; x<n_x-1; x++)
        {
            auto i_cell = x + y*n_x;
            neighbors[i_cell] = int_vec_t(3);
            neighbors[i_cell][0] = i_cell + n_x*plus_or_minus;
            neighbors[i_cell][1] = i_cell - 1; 
            neighbors[i_cell][2] = i_cell + 1;
        }
    };
    // Bottom-left corner
    auto connect_bounded_bl_corner = [&]()
    {
        auto i_cell = 0;
        neighbors[i_cell] = int_vec_t(2);
        neighbors[i_cell][0] = i_cell + 1; 
        neighbors[i_cell][1] = i_cell + n_x;
    };
    // Bottom-right corner
    auto connect_bounded_br_corner = [&]()
    {
        auto i_cell = n_x-1;
        neighbors[i_cell] = int_vec_t(2);
        neighbors[i_cell][0] = i_cell - 1;
        neighbors[i_cell][1] = i_cell + n_x;
    };
    // Top-left corner
    auto connect_bounded_tl_corner = [&]()
    {
        auto i_cell = (n_y-1)*n_x;
        neighbors[i_cell] = int_vec_t(2);
        neighbors[i_cell][0] = i_cell + 1; 
        neighbors[i_cell][1] = i_cell - n_x;
    };
    // Top-right corner
    auto connect_bounded_tr_corner = [&]()
    {
        auto i_cell = (n_x-1)+(n_y-1)*n_x;
        neighbors[i_cell] = int_vec_t(2);
        neighbors[i_cell][0] = i_cell - 1; 
        neighbors[i_cell][1] = i_cell - n_x;
    };

    /////////////////////////////////////////////

    // Step 1: Wire all the non-edge grid cells.
    connect_central_cells();

    // Step 2: Wire grid edge cells according to topology specs.
    switch (p.grid_topology)
    {
        case GridTopology::PERIODIC:
        {
            // Periodic grid topology in both x and y
            std::cout << "construct_2D_grid_multitopology: " 
                << "periodic" << std::endl;

            // Bottom row
            connect_periodic_x_edge_cells(0);
            // Top row
            connect_periodic_x_edge_cells(n_y-1);
            // Left column
            connect_periodic_y_edge_cells(0);
            // Right column
            connect_periodic_y_edge_cells(n_x-1);
            // Bottom-left corner
            connect_periodic_bl_corner();
            // Bottom-right corner
            connect_periodic_br_corner();
            // Top-left corner
            connect_periodic_tl_corner();
            // Top-right corner
            connect_periodic_tr_corner();

            return true;
        }
        case GridTopology::BOUNDED:
        {
            // Bounded grid topology in both x and y
            std::cout << "construct_2D_grid_multitopology: " 
                << "bounded" << std::endl;

            // Bottom row
            connect_bounded_x_edge_cells(0);
            // Top row
            connect_bounded_x_edge_cells(n_y-1);
            // Left column
            connect_bounded_y_edge_cells(0);
            // Right column
            connect_bounded_y_edge_cells(n_x-1);
            // Bottom-left corner
            connect_bounded_bl_corner();
            // Bottom-right corner
            connect_bounded_br_corner();
            // Top-left corner
            connect_bounded_tl_corner();
            // Top-right corner
            connect_bounded_tr_corner();

            return true;
        }
        default:
        {
            std::cout<< "construct_2D_grid_multitopology: " << "FAILED" << std::endl;
            return false;
        }
    }
}