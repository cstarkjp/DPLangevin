#
# Configuration file for Meson build system.
#

project(
    'DPLangevin',
    'cpp',
    version: '2025.10.11a1',
    default_options: ['cpp_std=c++11'],
)
add_project_arguments('-Ofast', '-g3', language : 'cpp')

cpp_sources = files(
    'src/langevin_construct_grid.cpp', 
    'src/langevin_grid_1d.cpp', 
    'src/langevin_grid_2d.cpp', 
    'src/langevin_ic.cpp', 
    'src/langevin_bc.cpp', 
    'src/langevin_set_coefficients.cpp', 
    'src/langevin_integrate_rungekutta.cpp', 
    'src/langevin_integrate_euler.cpp', 
    'src/langevin_utilities.cpp', 
    'src/application_dplangevin.cpp', 
    'src/application_dpsim_public.cpp', 
    'src/application_dpsim_private.cpp', 
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