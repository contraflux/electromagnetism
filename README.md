# electromagnetism
Simulation of electromagnetic waves written with Julia.

## Method
There are two simulations included in this package. The first uses the
Jefimenko-Feynman equations, and the second uses a finite time step method
and evolves Maxwell's equations using the Runge-Kutta 4 solver.

### Jefimenko-Feynman
To use this simulation, first run **jefimenko-feynman.jl** with the desired
parameters and point charge generating function at the top of the file. Then,
run **plot.jl**, varying colors, camera movement, and filename as needed.

### Causal-RK4
This file is a work in progress, but will run very similarly to the
Jefimenko-Feynman script once completed.

## Runs
A number of example simulations can be found in **/runs**

## Other
Three external scripts, **calculus.jl**, **physics.jl**, and
**runge_kutta.jl** contain helper functions for both the
Jefimenko-Feynman and Causal-RK4 simulations. The last script,
**pc_configs.jl** contains some interesting configurations
of point charges that can be inserted into either simulation.