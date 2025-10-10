/**
 * @file langevin_grid_2d.cpp
 * @brief Method for setting up a 2D grid for the model Langevin field.
 */

#include "general_core.hpp"

/**
* @details 
* Construct a connected 2D grid to be used for solving the evolution
* of a density field ρ(x,t).
* 
* Implements mixed grid edge topology: i.e., x and y grid edges 
* can separately be specified as "bounded" or "periodic".
* The "bounded" edge topology simply means those exterior grid cells are only
* connected to their neighboring grid-interior cells.
* The "periodic" edge topology means the grid cells along one edge are also 
* connected to those on the opposite edge. 
* If both x and y edges are periodic, the grid topology is toroidal.
* If only one edge is periodic, the grid topology is cylindrical.
* If both are bounded, the grid topology is a bounded plane.
*
* @param BaseLangevin integrator Parameters bundle.
*/
bool BaseLangevin::construct_2D_grid(const Parameters p)
{
    // Shorthand
    const auto n_x = p.n_x;
    const auto n_y = p.n_y;

    // Flattened grid vector each with a set of <=4 connection i_node elements.
    // Each i_node element will link to 1 of <=4 possible neighbor locations.
    // Along edges these sets will be reduced to 3 elements.
    // At corners these sets will be reduced to 2 elements.
    neighbors = std::vector<int_vec_t>(n_x*n_y, int_vec_t(0));

    auto i_from_xy = [&](int x, int y) -> int { return  x + y*n_x; };

    // Central cells
    auto wire_central_cell = [&](int x, int y)
    {
        // i_cell is the index of the flattened grid
        auto i_cell = i_from_xy(x, y);
        // Each cell has 4 neighbors.at(i_cell) indexes
        neighbors.at(i_cell).push_back(i_cell + n_x);  // Up:   i_cell + n_x // x + (y+1)*n_x;
        neighbors.at(i_cell).push_back(i_cell - n_x);  // Down: i_cell - n_x // x + (y-1)*n_x;
        neighbors.at(i_cell).push_back(i_cell + 1);    // Right: i_cell+1   (VMB: left)  // (x+1) + y*n_x;
        neighbors.at(i_cell).push_back(i_cell - 1);    // Left:  i_cell-1   (VMB: right) // (x-1) + y*n_x;

    };
    auto wire_central_cells = [&]()
    {
        for (auto y=1; y<n_y-1; y++)
        {
            for (auto x=1; x<n_x-1; x++)
            {
                wire_central_cell(x,y);
            }
        }
    };
    
    // Periodic
    auto wire_periodic_edge_cell_yplus = [&](int x, int y) 
    {
        auto i_cell = i_from_xy(x, y);
        auto i_yplus = (y < n_y-1) ? i_cell + n_x : x;
        neighbors.at(i_cell).push_back(i_yplus);   // Up
    };
    auto wire_periodic_edge_cell_yminus = [&](int x, int y) 
    {
        auto i_cell = i_from_xy(x, y);
        auto i_yminus = (y > 0) ? i_cell - n_x : x + (n_y-1)*n_x;
        neighbors.at(i_cell).push_back(i_yminus);  // Down
    };
    auto wire_periodic_edge_cell_xplus = [&](int x, int y) 
    {
        auto i_cell = i_from_xy(x, y);
        auto i_xplus = (x < n_x-1) ? i_cell + 1 : 0 + y*n_x;
        neighbors.at(i_cell).push_back(i_xplus);  // Right   (VMB: left)
    };
    auto wire_periodic_edge_cell_xminus = [&](int x, int y) 
    {
        auto i_cell = i_from_xy(x, y);
        auto i_xminus = (x > 0) ? i_cell - 1 : n_x-1 + y*n_x;
        neighbors.at(i_cell).push_back(i_xminus); // Left    (VMB: right)
    };
    auto wire_periodic_edge_cell = [&](int x, int y) 
    {
        wire_periodic_edge_cell_yplus(x, y);
        wire_periodic_edge_cell_yminus(x, y);
        wire_periodic_edge_cell_xplus(x, y);
        wire_periodic_edge_cell_xminus(x, y);        
    };
    auto wire_periodic_y_edges = [&](int x)
    {
        assert(x==0 or x==n_x-1);
        for (auto y=1; y<n_y-1; y++)
        {
            wire_periodic_edge_cell(x, y);
        }
    };
    auto wire_periodic_x_edges = [&](int y)
    {
        assert(y==0 or y==n_y-1);
        for (auto x=1; x<n_x-1; x++)
        {
            wire_periodic_edge_cell(x, y);
        }
    };

    // Bounded
    auto wire_bounded_y_edges = [&](int x)
    {
        assert(x==0 or x==n_x-1);
        auto plus_or_minus = (x==0) ? +1 : -1;
        for (auto y=1; y<n_y-1; y++)
        {
            auto i_cell = i_from_xy(x, y);
            neighbors.at(i_cell).push_back(i_cell + n_x*plus_or_minus);
            neighbors.at(i_cell).push_back(i_cell - n_x*plus_or_minus); 
            neighbors.at(i_cell).push_back(i_cell + plus_or_minus);
        }
    };
    auto wire_bounded_x_edges = [&](int y)
    {
        assert(y==0 or y==n_y-1);
        auto plus_or_minus = (y==0) ? +1 : -1;
        for (auto x=1; x<n_x-1; x++)
        {
            auto i_cell = i_from_xy(x, y);
            neighbors.at(i_cell).push_back(i_cell + n_x*plus_or_minus);
            neighbors.at(i_cell).push_back(i_cell - 1); 
            neighbors.at(i_cell).push_back(i_cell + 1);
        }
    };
    
    // Any topology
    auto wire_x_edges = [&](bool is_periodic_x_edge)
    {
        if (is_periodic_x_edge)
        {
            wire_periodic_x_edges(0);      // Bottom row
            wire_periodic_x_edges(n_y-1);  // Top row
        }
        else
        {
            wire_bounded_x_edges(0);        // Bottom row
            wire_bounded_x_edges(n_y-1);    // Top row
        }
    };
    auto wire_y_edges = [&](bool is_periodic_y_edge)
    {
        if (is_periodic_y_edge)
        {
            wire_periodic_y_edges(0);      // Left column
            wire_periodic_y_edges(n_x-1);  // Right column
        }
        else
        {
            wire_bounded_y_edges(0);        // Left column
            wire_bounded_y_edges(n_x-1);    // Right column
        }
    };

    // Corners
    auto wire_corners = [&](bool is_periodic_x_edge, bool is_periodic_y_edge)
    {
        int x, y, i_cell, i_yminus, i_xminus, i_yplus, i_xplus;

        // Bottom-left corner
        x = 0;
        y = 0;
        i_cell = i_from_xy(x, y);
        i_xplus = i_cell + 1;
        i_xminus = n_x-1 + y*n_x;
        i_yplus = i_cell + n_x;
        i_yminus = x + (n_y-1)*n_x;
        neighbors.at(i_cell).push_back(i_xplus);
        neighbors.at(i_cell).push_back(i_yplus);
        if (is_periodic_x_edge) { neighbors.at(i_cell).push_back(i_yminus); }
        if (is_periodic_y_edge) { neighbors.at(i_cell).push_back(i_xminus); }

        // Top-left corner
        x = 0;
        y = n_y-1;
        i_cell = i_from_xy(x, y);
        i_xplus = i_cell + 1;
        i_xminus = n_x-1 + y*n_x;
        i_yplus = x;
        i_yminus = i_cell - n_x;
        neighbors.at(i_cell).push_back(i_xplus);
        neighbors.at(i_cell).push_back(i_yminus);
        if (is_periodic_x_edge) { neighbors.at(i_cell).push_back(i_yplus); }
        if (is_periodic_y_edge) { neighbors.at(i_cell).push_back(i_xminus); }

        // Bottom-right corner
        x = n_x-1;
        y = 0;
        i_cell = i_from_xy(x, y);
        i_xplus = 0 + y*n_x;
        i_xminus = i_cell - 1;
        i_yplus = i_cell + n_x;
        i_yminus = x + (n_y-1)*n_x;
        neighbors.at(i_cell).push_back(i_xminus);
        neighbors.at(i_cell).push_back(i_yplus);
        if (is_periodic_x_edge) { neighbors.at(i_cell).push_back(i_yminus); }
        if (is_periodic_y_edge) { neighbors.at(i_cell).push_back(i_xplus); }

        // Top-right corner
        x = n_x-1;
        y = n_y-1;
        i_cell = i_from_xy(x, y);
        i_xplus = 0 + y*n_x;
        i_xminus = i_cell - 1;
        i_yplus = x;
        i_yminus = i_cell - n_x;
        neighbors.at(i_cell).push_back(i_xminus);
        neighbors.at(i_cell).push_back(i_yminus);
        if (is_periodic_x_edge) { neighbors.at(i_cell).push_back(i_yplus); }
        if (is_periodic_y_edge) { neighbors.at(i_cell).push_back(i_xplus); }
    };
    
    // Whole grid
    auto wire_grid = [&](gt_vec_t grid_topologies)
    {
        auto is_periodic_x_edge = (grid_topologies[0]==GridTopology::PERIODIC);
        auto is_periodic_y_edge = (grid_topologies[1]==GridTopology::PERIODIC);
        wire_x_edges(is_periodic_x_edge);
        wire_y_edges(is_periodic_y_edge);
        wire_corners(is_periodic_x_edge, is_periodic_y_edge);
    };

    /////////////////////////////////////////////

    // Step 1: Wire all the non-edge grid cells.
    wire_central_cells();

    // Step 2: Wire grid edge cells according to topology specs.
    wire_grid(p.grid_topologies);

    return true;
}