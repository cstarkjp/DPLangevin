#!/usr/bin/env python

import numpy as np
from numpy.typing import NDArray
import dplvn # type: ignore

print()
result: NDArray = dplvn.dp(
    n_cells=64, #4096 
    linear=1.0, quadratic=2.0, diffusion=0.1, noise=1.0, 
    # t_max=100.0, dx=0.5, dt=0.01,
    # grid_dimension=dplvn.D2,
    # initial_condition=dplvn.RANDOM_UNIFORM,
    # boundary_condition=dplvn.PERIODIC,
)
print()
print( type(result) )
print( result.shape )
print( result[-10:] )
print()

epochs: NDArray = result[:, 0]
mean_densities: NDArray = result[:, 1]
print()
print( epochs[-10:] ) 
print( mean_densities[-10:] )