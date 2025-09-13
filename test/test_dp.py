#!/usr/bin/env python

import numpy as np
from numpy.typing import NDArray
import dplvn # type: ignore

print()
print(f"dplvn version:  {dplvn.__version__}")
print()
sim = dplvn.dp(
    linear=1.0, quadratic=2.0, diffusion=0.1, noise=1.0, 
    t_max=100.0, 
    # dx=0.5, dt=0.01
    # grid_dimension=dplvn.D1,
    # grid_size=(4096,),
    grid_dimension=dplvn.D2,
    grid_size=(8,8,),
    grid_topology= dplvn.PERIODIC,
    boundary_condition=dplvn.FLOATING,
    integration_method=dplvn.RUNGE_KUTTA
)
print()
print()

# print(sim.get())
print(sim.get_epochs())
print()
print(sim.get_mean_densities())
print()

# print()
# print( type(result) )
# print( result.shape )
# print( result[-10:] )
# print()

# epochs: NDArray = result[:, 0]
# mean_densities: NDArray = result[:, 1]
# print()
# print( epochs[-10:] ) 
# print( mean_densities[-10:] )

# arrays = dplvn.assign_results()
# print(arrays.get())
# arrays.setName()
# print(arrays.getName())