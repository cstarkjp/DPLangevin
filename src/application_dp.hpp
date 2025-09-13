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


class Results {
private:
    std::string m_name;
    dbl_vec_t epochs;
    dbl_vec_t mean_densities;
    results_t results;
public:
    Results(
        const std::string& name
    )
    {
        m_name = name;
        epochs = {0.0, 1, 2, 3};
        mean_densities = {10,20,30,40.0};
    }
    results_t get() const { return results; }

    void set_arrays() 
    {
        epochs = {0.0, 1, 2, 3};
        mean_densities = {10,20,30,40.0};
    }

    void prepare_return_array(
        // const dbl_vec_t& cell_density,
        const dbl_vec_t& epochs, const dbl_vec_t& mean_densities
    )
    {
        results_t results_array({static_cast<int>(epochs.size()), 2});
        auto array_proxy = results_array.mutable_unchecked();
        for (auto i=0; i<epochs.size(); i++)
        {
            array_proxy(i, 0) = epochs[i];
            array_proxy(i, 1) = mean_densities[i];
        };
        results = results_array;
    }
};

// results_t
auto dp(
    const double linear, const double quadratic, 
    const double diffusion, const double noise, 
    const double t_max, const double dx, const double dt, const int random_seed,
    const GridDimension grid_dimension, 
    const int_vec_t& grid_size,
    const GridTopology grid_topology,
    const BoundaryCondition boundary_condition,
    const InitialCondition initial_condition,
    const IntegrationMethod integration_method
) -> Results;

#endif
