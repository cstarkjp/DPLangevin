project(
    'DPLangevin',
    'cpp',
    version: '2025.09.14a2',
    default_options: [
        'cpp_std=c++11',
    ],
)

add_project_arguments('-Ofast', '-g3', language : 'cpp')
py = import('python').find_installation(pure: false)
pybind11_dep = dependency('pybind11')

cpp_sources = files(
    'src/langevin_grid_1d.cpp', 
    'src/langevin_grid_2d.cpp', 
    'src/langevin_initial_condition.cpp', 
    'src/langevin_initialization.cpp', 
    'src/langevin_integration_rungekutta.cpp', 
    'src/langevin_integration_euler.cpp', 
    'src/langevin_utils.cpp', 
    'src/application_dplangevin.cpp', 
    'src/application_dpsim_public.cpp', 
    'src/application_dpsim_private.cpp', 
    'src/applications_pybind.cpp'
)

py.extension_module(
    'dplvn',
    cpp_sources,
    install: true,
    dependencies : [pybind11_dep]
)