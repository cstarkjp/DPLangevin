#!/usr/bin/env python

import numpy as np
from numpy.typing import NDArray
import dplvn # type: ignore

print()
print(f"dplvn version:  {dplvn.__version__}")
print()

sim_dp = dplvn.SimDP(
    linear=1.0, quadratic=2.0, 
    diffusion=0.1, noise=1.0, 
    t_max=50.0-1e-10, 
    # dx=0.5, dt=0.01
    # grid_dimension=dplvn.D1,
    # grid_size=(4096,),
    grid_dimension=dplvn.D2,
    grid_size=(8,8,),
    # grid_size=(64,64,),
    grid_topology= dplvn.BOUNDED,
    boundary_condition=dplvn.FLOATING,
    initial_condition=dplvn.RANDOM_UNIFORM,
    integration_method=dplvn.RUNGE_KUTTA
)
print(dir(sim_dp))
if not sim_dp.initialize():
    raise Exception("Failed to initialize sim")
if not sim_dp.run():
    raise Exception("Failed to run sim")
print(sim_dp.get_epochs())
print(sim_dp.get_mean_densities())
print(sim_dp.get_density())

