#!/usr/bin/env python

import numpy as np
from numpy.typing import NDArray
import dplvn # type: ignore

result = dplvn.dp(
    n_cells=64, #4096 
    linear=1.0, quadratic=2.0, diffusion=0.1, noise=1.0, 
    # t_max=100.0, dx=0.5, dt=0.01
)
print( type(result) )
print( result.shape )
print( result[-10:] )

epochs: NDArray = result[:, 0]
mean_densities: NDArray = result[:, 1]
print( epochs[-10:], mean_densities[-10:] )