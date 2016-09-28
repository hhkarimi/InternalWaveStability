# PDEevolution
Solves the partial differential equations that describe the time-evolution of the underlying internal gravity wave beam subject to subharmonic disturbances outline in chapter 3 of [my PhD thesis](https://dspace.mit.edu/handle/1721.1/100060).

## `solvePDE.m`
Main function that sends the PDE defined in `mlinesPDE.m` through an adaptive Runge-Kutta intergration algorithm to solve for the wave profiles.  The solution is then plotted in various forms of visualizations.  A mesh/surface plot presents the physical phenomena of instability most clearly.

## `mlinesPDE.m`
Partial differential equation describing the evolution of the underlying internal gravity wave beam with a pair of subharmonic disturbances.  The resulting nonlinear system, due to triad interations, is discretized in space via finite differences.  At each grid point, an ordinary differential equation in time is written and passed into an adaptive RK45 integration scheme.  This approach is known as the "method of lines".
