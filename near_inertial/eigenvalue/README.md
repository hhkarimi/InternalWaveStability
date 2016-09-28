# eigenvalue
Solves the eigenvalue problem that determines the instability of internal gravtiy wave beams to PSI.  This is an implementation of the mathematical analysis performed in chapter 3 of [my PhD thesis](https://dspace.mit.edu/handle/1721.1/100060).

## `main.m`
Loads an initial guess for eigenvalues (usually found from a finite difference scheme) to load into optimization problem solved by Newton's method.  The associated boundary value problem is solved using the shooting method (guessing boundary values at one end and minimizing the error at the other by integrating with standard initial value tools).  `main.m` calls the optimization routine and outputs the eigenvalues as they depend on the parameter kappa.

## `beam_profile.m`
Defines the nondimensional, cross-sectional profile of the beam.  Typically taken to be a Gaussian or a tanh-tanh (smoothened top-hat) function.

## `ddn.m`
