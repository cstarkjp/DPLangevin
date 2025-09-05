project(
    'DPLangevin',
    'cpp',
    version: '2025.9.2',
    default_options: [
        'cpp_std=c++11',
    ],
)

# if get_option('buildtype') == 'custom'
add_project_arguments('-Ofast', '-g3', language : 'cpp')
# endif


py = import('python').find_installation(pure: false)
pybind11_dep = dependency('pybind11')
# ninja_dep = dependency('ninja')
cpp_sources = files(
    'src/integration.cpp',
    'src/initialization.cpp', 
    'src/grid1d.cpp', 
    'src/grid2d.cpp', 
    'src/network.cpp', 
    'src/ic.cpp', 
    'src/rungekutta.cpp', 
    'src/eulermethod.cpp', 
    'src/utils.cpp', 
    'src/nonlinear.cpp', 
    'src/demo.cpp', 
    'src/wrapper.cpp'
)
module_name = 'dplvn'

py.extension_module(
    module_name,
    cpp_sources,
    install: true,
    dependencies : [pybind11_dep],
)