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
typedef py::array_t<double, py::array::c_style> py_array_t;


class Results {
private:
    std::string m_name;
    py_array_t results;
    py_array_t return_epochs;
    py_array_t return_mean_densities;
    int n_epochs;
public:
    Results(
        const std::string& name
    )
    {
        m_name = name;
    }

    void prepare_return_arrays(
        // const dbl_vec_t& cell_density,
        const dbl_vec_t& epochs, const dbl_vec_t& mean_densities
    )
    {
        n_epochs = static_cast<int>(epochs.size());

        py_array_t results_array({n_epochs, 2});
        auto array_proxy = results_array.mutable_unchecked();
        for (auto i=0; i<n_epochs; i++)
        {
            array_proxy(i, 0) = epochs[i];
            array_proxy(i, 1) = mean_densities[i];
        };
        results = results_array;

        py_array_t epochs_array(n_epochs);
        auto epochs_proxy = epochs_array.mutable_unchecked();
        for (auto i=0; i<n_epochs; i++)
        {
            epochs_proxy(i) = epochs[i];
        };
        return_epochs = epochs_array;

        py_array_t mean_densities_array(n_epochs);
        auto mean_densities_proxy = mean_densities_array.mutable_unchecked();
        for (auto i=0; i<n_epochs; i++)
        {
            mean_densities_proxy(i) = mean_densities[i];
        };
        return_mean_densities = mean_densities_array;
    }

    py_array_t get() const { return results; }
    py_array_t get_epochs() const { return return_epochs; }
    py_array_t get_mean_densities() const { return return_mean_densities; }
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
