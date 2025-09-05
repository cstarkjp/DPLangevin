// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#include <pybind11/numpy.h>
#include "declarations.hpp"
#include "application.hpp"

PYBIND11_MODULE(dplvn, m) {
    // The main function
    m.doc() = "'Dornic' operator-splitting method of integrating DP-type \
Langevin equations"; 
    m.def(
        "demo", 
        &demo,
        py::arg("linear") = 1.0, 
        py::arg("quadratic") = 2.0, 
        py::arg("diffusion") = 0.1,
        py::arg("noise") = 1.0,
        py::arg("n_cells") = 4096, 
        py::arg("t_max") = 100.0,
        py::arg("dx") = 0.5,
        py::arg("dt") = 0.01,
        py::arg("random_seed") = 1, //45345456, 
        "Demo application of the Dornic method"
    );
}