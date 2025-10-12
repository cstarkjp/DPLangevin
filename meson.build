#
# Configuration file for Meson build system.
#

project(
    'DPLangevin',
    'cpp',
    version: '2025.10.12a1',
    default_options: ['cpp_std=c++11'],
)
add_project_arguments('-Ofast', '-g3', language : 'cpp')

cpp_sources = files(
    'src/langevin_construct_grid.cpp', 
    'src/langevin_construct_grid1d.cpp', 
    'src/langevin_construct_grid2d.cpp', 
    'src/langevin_ic.cpp', 
    'src/langevin_bc.cpp', 
    'src/langevin_prepare.cpp', 
    'src/langevin_integrate_rungekutta.cpp', 
    'src/langevin_integrate_euler.cpp', 
    'src/langevin_utilities.cpp', 
    'src/dplangevin.cpp', 
    'src/sim_dplangevin.cpp', 
    'src/sim_dplangevin_private.cpp', 
    'src/sim_dplangevin_utilities.cpp',
    'src/wrapper_pybind.cpp'
)

pybind11_dep = dependency('pybind11')
py = import('python').find_installation(pure: false)
py.extension_module(
    'dplvn',
    cpp_sources,
    install: true,
    dependencies : [pybind11_dep]
)