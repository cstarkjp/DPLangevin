project(
    'DPLangevin',
    'cpp',
    version: '2025.09.8a9',
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

# py_sources = files(
#     'src/__init__.py',
#     'src/initialize.py'
# )
# py.install_sources(
#     py_sources,
#     subdir: 'dplvn'
# )

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

# https://mesonbuild.com/meson-python/how-to-guides/shared-libraries.html
# example_lib = shared_library(
#     'dplvn',
#     cpp_sources,
#     install: true,
#     install_dir: py.get_install_dir() / 'mypkg/subdir',
# )

py.extension_module(
    'dplvn',
    cpp_sources,
    install: true,
    # subdir: 'dplvn',
    # install_dir: py.get_install_dir(),
    dependencies : [pybind11_dep]
)