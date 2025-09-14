# rm -rf build; meson setup build; meson compile -C build

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
    'src/integration_rungekutta.cpp', 
    'src/integration_euler.cpp', 
    'src/grid_1d.cpp', 
    'src/grid_2d.cpp', 
    'src/initial_condition.cpp', 
    'src/initialization.cpp', 
    'src/utils.cpp', 
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