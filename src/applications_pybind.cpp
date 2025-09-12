// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 
// 
#include <pybind11/numpy.h>
#include "core.hpp"
#include "application_dp.hpp"
// Essential for STL container conversions
#include <pybind11/stl.h> 

PYBIND11_MODULE(dplvn, module)
{
    module.attr("__version__") = "2025.09.12a5";
    module.doc() = 
        "'Dornic' operator-splitting method of integrating DP-type Langevin equations"; 

    py::enum_<GridDimension>(module, "GridDimension")
        .value("D1", GridDimension::D1)
        .value("D2", GridDimension::D2)
        .value("D3", GridDimension::D3)
        .value("D4", GridDimension::D4)
        .export_values();

    py::enum_<GridTopology>(module, "GridTopology")
        .value("BOUNDED", GridTopology::BOUNDED)
        .value("PERIODIC", GridTopology::PERIODIC)
        .export_values();

    py::enum_<BoundaryCondition>(module, "BoundaryCondition")
        .value("FLOATING", BoundaryCondition::FLOATING)
        .value("FIXED_VALUE", BoundaryCondition::FIXED_VALUE)
        .value("FIXED_FLUX", BoundaryCondition::FIXED_FLUX)
        .export_values();

    py::enum_<InitialCondition>(module, "InitialCondition")
        .value("RANDOM_UNIFORM", InitialCondition::RANDOM_UNIFORM)
        .value("RANDOM_GAUSSIAN", InitialCondition::RANDOM_GAUSSIAN)
        .value("CONSTANT_VALUE", InitialCondition::CONSTANT_VALUE)
        .value("SINGLE_SEED", InitialCondition::SINGLE_SEED)
        .export_values();

    py::enum_<IntegrationMethod>(module, "IntegrationMethod")
        .value("EULER", IntegrationMethod::EULER)
        .value("RUNGE_KUTTA", IntegrationMethod::RUNGE_KUTTA)
        .export_values();
        
    module.def(
        "dp", 
        &dp,
        py::arg("linear") = 1.0, 
        py::arg("quadratic") = 2.0, 
        py::arg("diffusion") = 0.1,
        py::arg("noise") = 1.0,
        py::arg("t_max") = 100.0,
        py::arg("dx") = 0.5,
        py::arg("dt") = 0.01,
        py::arg("random_seed") = 1,
        py::arg("grid_dimension") = GridDimension::D2,
        py::arg("grid_size") = py::none(),
        py::arg("grid_topology") = GridTopology::BOUNDED,
        py::arg("boundary_condition") = BoundaryCondition::FLOATING,
        py::arg("initial_condition") = InitialCondition::RANDOM_UNIFORM,
        py::arg("integration_method") = IntegrationMethod::RUNGE_KUTTA,
        "Demo application of the Dornic method"
    );
}