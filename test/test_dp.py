#!/usr/bin/env python

import numpy as np
from numpy.typing import NDArray
import dplvn # type: ignore

bold = lambda str: ("\033[1m" + str + "\033[0m")

print()
print(bold(f"dplvn version:  {dplvn.__version__}"))
print()

Δt: float = 0.01
sim = dplvn.SimDP(
    linear=1.0, quadratic=2.0, 
    diffusion=0.1, noise=1.0, 
    t_final=50.0-1e-10, 
    # dx=0.5, 
    dt=Δt,
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

if not sim.initialize():
    raise Exception("Failed to initialize sim")

n_segments: int = 5
n_epochs: int = sim.get_n_epochs()
n_segment_epochs: int = (n_epochs-1) // n_segments
if (n_segment_epochs*n_segments+1)!=n_epochs:
    raise Exception(
        f"Failed to segment sim with {n_epochs} epochs "
        + "into {n_segments} segment(s)"
    )

print(bold(f"Integrating:  {n_epochs} epochs in {n_segments} segment(s)"))
print()
# Clean this repetition up!
i_segment: int = 0
if not sim.process():
    raise Exception("Failed to process sim results")
i_epoch: int = sim.get_i_epoch()
t_epoch: float = sim.get_t_epoch()
print(bold(f"segment={i_segment}/{n_segments}  i={i_epoch-1}"))
print(f"t epochs:  {sim.get_t_epochs()}")
print(f"mean densities:  {sim.get_mean_densities()}")
print("cell density grid:")
print(np.round(sim.get_density().T, 2))
print()
for i_segment in range(n_segments):
    if not sim.run(n_segment_epochs):
        raise Exception("Failed to run sim")
    if not sim.process():
        raise Exception("Failed to process sim results")
    i_epoch = sim.get_i_epoch()
    t_epoch = sim.get_t_epoch()
    print(bold(
        f"segment={i_segment+1}/{n_segments}  "
        + f"i={i_epoch-1} t={np.round(t_epoch-Δt,5)}"
    ))
    print(f"t epochs:  {sim.get_t_epochs()}")
    print(f"mean densities:  {sim.get_mean_densities()}")
    print("cell density grid:")
    print(np.round(sim.get_density().T, 2))
    print()
