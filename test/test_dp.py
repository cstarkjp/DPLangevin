#!/usr/bin/env python

import numpy as np
from numpy.typing import NDArray
import dplvn # type: ignore

bold = lambda str: ("\033[1m" + str + "\033[0m")

print()
print(bold(f"dplvn version:  {dplvn.__version__}"))
print()

sim_dp = dplvn.SimDP(
    linear=1.0, quadratic=2.0, 
    diffusion=0.1, noise=1.0, 
    t_final=50.0-1e-10, 
    # dx=0.5, dt=0.01,
    random_seed=1,
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
print()

if not sim_dp.initialize():
    raise Exception("Failed to initialize sim")

n_segments: int = 5
n_epochs: int = sim_dp.get_n_epochs()
n_segment_epochs: int = n_epochs // n_segments
if (n_segment_epochs*n_segments+1)!=n_epochs:
    raise Exception(
        f"Failed to segment sim with {n_epochs} epochs into {n_segments}"
    )

print(bold("Integrating:"))
print()
for i_segment in range(n_segments):
    print(bold(f"i={i_segment+1}/{n_segments}"))
    if not sim_dp.run(n_segment_epochs):
        raise Exception("Failed to run sim")
    if not sim_dp.finalize():
        raise Exception("Failed to finalize sim")
    print(f"epochs:  {sim_dp.get_epochs()}")
    print(f"mean_densities:  {sim_dp.get_mean_densities()}")
    print("cell density grid (at final t):")
    print(np.round(sim_dp.get_density().T, 2))
    print()
print()