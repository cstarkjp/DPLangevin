// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef APPLICATION_DP_HPP
#define APPLICATION_DP_HPP

namespace py = pybind11;
typedef py::array_t<double, py::array::c_style> results_t;

results_t dp(
    const double linear, const double quadratic, 
    const double diffusion, const double noise, 
    const double t_max, const double dx, const double dt, const int random_seed,
    const GridDimension grid_dimension, 
    const int_vec_t& grid_size,
    const GridTopology grid_topology,
    const BoundaryCondition boundary_condition,
    const InitialCondition initial_condition,
    const IntegrationMethod integration_method
);

class Results {
private:
    std::string m_name;
    dbl_vec_t test;
    // dbl_vec_t mean_densities(10, 0.0);
public:
    Results(
        const std::string& name)
    {
        m_name = name;
        test = {1, 2, 3};
    }
    std::string getName() const { return m_name; }
    void setName(const std::string& name) { m_name = name; }
};


#endif