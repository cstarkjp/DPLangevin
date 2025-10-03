/**
 * @file wrapper_pybind.cpp
 * @brief Pybind11 wrapper between C++ and Python for DP Langevin application.
 */

#include <pybind11/numpy.h>
// Essential for STL container conversions
#include <pybind11/stl.h> 
#include "general_core.hpp"
#include "application_dplangevin.hpp"
#include "application_dpsim.hpp"

PYBIND11_MODULE(dplvn, module)
{
    module.attr("__version__") = "2025.10.03c";
    module.doc() = 
        "Operator-splitting method of integrating DP-type Langevin equations"; 

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
        
    py::class_<SimDP>(module, "SimDP")
        .def(
            py::init<
                double, double, 
                double, double, 
                double, double, double,
                int, 
                GridDimension,
                int_vec_t&,
                GridTopology,
                BoundaryCondition,
                InitialCondition,
                IntegrationMethod
            >(),
            "Simulation of DP Langevin equation",
            py::arg("linear") = 1.0, 
            py::arg("quadratic") = 2.0, 
            py::arg("diffusion") = 0.1,
            py::arg("noise") = 1.0,
            py::arg("t_final") = 100.0,
            py::arg("dx") = 0.5,
            py::arg("dt") = 0.01,
            py::arg("random_seed") = 1,
            py::arg("grid_dimension") = GridDimension::D2,
            py::arg("grid_size") = int_vec_t(4),
            py::arg("grid_topology") = GridTopology::BOUNDED,
            py::arg("boundary_condition") = BoundaryCondition::FLOATING,
            py::arg("initial_condition") = InitialCondition::RANDOM_UNIFORM,
            py::arg("integration_method") = IntegrationMethod::RUNGE_KUTTA
        )
        .def("initialize", &SimDP::initialize)
        .def("run", &SimDP::run)
        .def("process", &SimDP::process)
        .def("get_n_epochs", &SimDP::get_n_epochs)
        .def("get_i_next_epoch", &SimDP::get_i_next_epoch)
        .def("get_i_current_epoch", &SimDP::get_i_current_epoch)
        .def("get_t_next_epoch", &SimDP::get_t_next_epoch)
        .def("get_t_current_epoch", &SimDP::get_t_current_epoch)
        .def("get_t_epochs", &SimDP::get_t_epochs)
        .def("get_mean_densities", &SimDP::get_mean_densities)
        .def("get_density", &SimDP::get_density);

}