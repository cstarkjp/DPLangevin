/**
 * @file langevin_grid_2d_mixedtopology.cpp
 * @brief Method for setting up a 2D grid for the model Langevin field.
 */

#include "general_core.hpp"

//! Utility to allow switch-case with GridTopology tuple (overloaded)
constexpr unsigned long pack(const gt_vec_t& gt) {
    return (
        (static_cast<unsigned long>(gt[0]) << 8) | 
        static_cast<unsigned long>(gt[1])
    );
}
//! Utility to allow switch-case with GridTopology tuple (overloaded)
constexpr unsigned long pack(GridTopology a, GridTopology b) {
    return (
        (static_cast<unsigned long>(a) << 8) | 
        static_cast<unsigned long>(b)
    );
}

/**
* @details Construct 2D density field ρ(x,t) grid and corresponding cell-cell topologies
* 
* Implements mixed grid edge topology: i.e., x and y grid edges 
* can separately be specified as "bounded" or "periodic".
*
* The "bounded" edge topology simply means those exterior grid cells are only
* connected to their neighboring grid-interior cells.
*
* The "periodic" edge topology means the grid cells along one edge are also 
* connected to those on the opposite edge. 
*
* If both x and y edges are periodic, the grid topology is toroidal.
* If only one edge is periodic, the grid topology is cylindrical.
* If both are bounded, the grid topology is a bounded plane.
*
* @param BaseLangevin integrator Parameters bundle.
*/
bool BaseLangevin::construct_2D_grid_mixedtopology(const Parameters p)
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
    auto wire_central_cell = [&](int x, int y)
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
        auto i_edge_cell = x + y*n_x;
        auto i_yplus = (y < n_y-1) ? i_edge_cell + n_x : x;
        neighbors[i_edge_cell][0] = i_yplus;   // Up
    };
    auto wire_periodic_edge_cell_yminus = [&](int x, int y) 
    {
        auto i_edge_cell = x + y*n_x;
        auto i_yminus = (y > 0) ? i_edge_cell - n_x : x + (n_y-1)*n_x;
        neighbors[i_edge_cell][1] = i_yminus;  // Down
    };
    auto wire_periodic_edge_cell_xplus = [&](int x, int y) 
    {
        auto i_edge_cell = x + y*n_x;
        auto i_xplus = (x < n_x-1) ? i_edge_cell + 1 : 0 + y*n_x;
        neighbors[i_edge_cell][2] = i_xplus;  // Right   (VMB: left)
    };
    auto wire_periodic_edge_cell_xminus = [&](int x, int y) 
    {
        auto i_edge_cell = x + y*n_x;
        auto i_xminus = (x > 0) ? i_edge_cell - 1 : n_x-1 + y*n_x;
        neighbors[i_edge_cell][3] = i_xminus; // Left    (VMB: right)
    };
    auto wire_periodic_edge_cells = [&](int x, int y) 
    {
        wire_periodic_edge_cell_yplus(x, y);
        wire_periodic_edge_cell_yminus(x, y);
        wire_periodic_edge_cell_xplus(x, y);
        wire_periodic_edge_cell_xminus(x, y);        
    };
    // Left and right edges, loop over y cells
    auto wire_periodic_y_edge_cells = [&](int x)
    {
        assert(x==0 or x==n_x-1);
        for (auto y=1; y<n_y-1; y++)
        {
            wire_periodic_edge_cells(x, y);
        }
    };
    // Bottom and top edges, loop over x cells
    auto wire_periodic_x_edge_cells = [&](int y)
    {
        assert(y==0 or y==n_y-1);
        for (auto x=1; x<n_x-1; x++)
        {
            wire_periodic_edge_cells(x, y);
        }
    };
    // Corners
    auto wire_periodic_corner = [&](int x, int y)
    {
        wire_periodic_edge_cell_yplus(x, y);
        wire_periodic_edge_cell_yminus(x, y);
        wire_periodic_edge_cell_xplus(x, y);
        wire_periodic_edge_cell_xminus(x, y);
    };

    // Bounded
    // Left and right edges, loop over y cells
    auto wire_bounded_y_edge_cells = [&](int x)
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
    auto wire_bounded_x_edge_cells = [&](int y)
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
    // Corners
    auto wire_bounded_corner = [&](int x, int y)
    {
        auto x_plus_or_minus = (x==0) ? +1 : -1;
        auto y_plus_or_minus = (y==0) ? +1 : -1;
        auto i_cell = x + y*n_x;
        neighbors[i_cell] = int_vec_t(2);
        neighbors[i_cell][0] = i_cell + 1*x_plus_or_minus; 
        neighbors[i_cell][1] = i_cell + n_x*y_plus_or_minus;
    };

    /////////////////////////////////////////////

    // Step 1: Wire all the non-edge grid cells.
    wire_central_cells();

    // Step 2: Wire grid edge cells according to topology specs.
    auto grid_topologies = pack(p.grid_topologies);
    switch (grid_topologies) 
    {
        case pack(GridTopology::PERIODIC, GridTopology::PERIODIC):
        {
            // Periodic grid topology in both x and y
            std::cout 
                << "construct_2D_grid_mixedtopology: " 
                << "x:periodic, y:periodic"
                << std::endl;

            wire_periodic_x_edge_cells(0);      // Bottom row
            wire_periodic_x_edge_cells(n_y-1);  // Top row
            wire_periodic_y_edge_cells(0);      // Left column
            wire_periodic_y_edge_cells(n_x-1);  // Right column
            wire_periodic_corner(0, 0);         // Bottom-left corner
            wire_periodic_corner(n_x-1, 0);     // Bottom-right corner
            wire_periodic_corner(0, n_y-1);     // Top-left corner
            wire_periodic_corner(n_x-1, n_y-1); // Top-right corner

            return true;
        }
        case pack(GridTopology::BOUNDED, GridTopology::BOUNDED):
        {
            // Bounded grid topology in both x and y
            std::cout 
                << "construct_2D_grid_mixedtopology: " 
                << "x:bounded, y:bounded"
                << std::endl;

            wire_bounded_x_edge_cells(0);        // Bottom row
            wire_bounded_x_edge_cells(n_y-1);    // Top row
            wire_bounded_y_edge_cells(0);        // Left column
            wire_bounded_y_edge_cells(n_x-1);    // Right column
            wire_bounded_corner(0, 0);           // Bottom-left corner
            wire_bounded_corner(n_x-1, 0);       // Bottom-right corner
            wire_bounded_corner(0, n_y-1);       // Top-left corner
            wire_bounded_corner(n_x-1, n_y-1);   // Top-right corner

            return true;
        }
        case pack(GridTopology::BOUNDED, GridTopology::PERIODIC):
        {
            // Periodic in x direction
            // Bounded along x edges
            // Periodic along y edges
            std::cout 
                << "construct_2D_grid_mixedtopology: " 
                << "x:bounded, y:periodic"
                << std::endl;

            wire_bounded_x_edge_cells(0);        // Bottom row
            wire_bounded_x_edge_cells(n_y-1);    // Top row
            wire_periodic_y_edge_cells(0);       // Left column
            wire_periodic_y_edge_cells(n_x-1);   // Right column
            // TBD: rewire for bounded directions
            wire_periodic_corner(0, 0);          // Bottom-left corner
            wire_periodic_corner(n_x-1, 0);      // Bottom-right corner
            wire_periodic_corner(0, n_y-1);      // Top-left corner
            wire_periodic_corner(n_x-1, n_y-1);  // Top-right corner

            return true;
        }
        case pack(GridTopology::PERIODIC, GridTopology::BOUNDED):
        {
            // Periodic in y direction
            // Periodic along x edges
            // Bounded along y edges
            std::cout 
                << "construct_2D_grid_mixedtopology: " 
                << "x:periodic, y:bounded"
                << std::endl;

            wire_periodic_x_edge_cells(0);       // Bottom row
            wire_periodic_x_edge_cells(n_y-1);   // Top row
            wire_bounded_y_edge_cells(0);        // Left column
            wire_bounded_y_edge_cells(n_x-1);    // Right column
            // TBD: rewire for bounded directions
            wire_periodic_corner(0, 0);          // Bottom-left corner
            wire_periodic_corner(n_x-1, 0);      // Bottom-right corner
            wire_periodic_corner(0, n_y-1);      // Top-left corner
            wire_periodic_corner(n_x-1, n_y-1);  // Top-right corner

            return true;
        }
        default:
        {
            std::cout
                << "construct_2D_grid_mixedtopology: " 
                << "FAILED "
                << std::hex << grid_topologies
                << std::endl;
            return false;
        }
    }
}