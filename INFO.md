
# Program design

The structure of the DP/APT Langevin-equation integrator package is broadly as follows 
(access the pertinent source `C++` files 
[here](https://github.com/cstarkjp/DPLangevin/tree/main/src/)).

At the top level, there is a wrapper file called [`wrapper_pybind.cpp`](https://github.com/cstarkjp/DPLangevin/tree/main/src/wrapper_pybind.cpp) that uses `pybind11` to link the `C++` code to a Python runtime. It `#include`s the "application" and core header files.

At the lower level, the code is split into three groups, each denoted by one of three file prefixes: (1) `application_`, (2) `langevin_` or (3) `general_`:

   1. The `application_` source files implement specific uses of the operator-splitting integration method. For now, the only implemented application is to solve the directed-percolation (DP) Langevin equation. 
   This implementation is split into a "simulation" part (`application_dpsim_` files) and an "integrator" part (`application_dplangevin_` files). 

       The `application_dpsim_` files provide a `SimDP` class, made available through the wrapper at the Python level, required to manage and execute DP Langevin model integration.  Each instance of the `SimDP` class instantiates a `DPLangevin` class integrator to do the hard work of numerical integration of the stochastic differential equation.

       The `application_dplangevin_` files define this `DPLangevin` integrator class. They inherit the general `Langevin` integrator class and implement several methods left undefined by that parent; most important, they define methods implementing the particular functional form of the directed-percolation Langevin equation and its corresponding nonlinear, deterministic integration step in the split operator scheme.


   2. The `langevin_` source files provide the base `Langevin` class that implements the operator-splitting integration method in a fairly general fashion. Grid geometry and topology, boundary conditions, initial conditions, the integration scheme, and a general form of the Langevin equation are all coded here. Some of these methods are heavily altered versions of the Villa-Martín and Buendían code; others remain very similar to their original implementations.

   3. The `general_` source files provide the basic stuff used throughout the code, notably the `typedef`s, `struct`s, `enum`s, and macros.