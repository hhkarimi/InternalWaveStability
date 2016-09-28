# eigenvalue
Solves the eigenvalue problem that determines the instability of internal gravtiy wave beams to PSI.  This is an implementation of the mathematical analysis performed in chapter 3 of [my PhD thesis](https://dspace.mit.edu/handle/1721.1/100060).

## `main.m`
Loads an initial guess for eigenvalues (usually found from a finite difference scheme) to load into optimization problem solved by Newton's method.  The associated boundary value problem is solved using the shooting method (guessing boundary values at one end and minimizing the error at the other by integrating with standard initial value tools).  `main.m` calls the optimization routine and outputs the eigenvalues as they depend on the parameter kappa.

## `beam_profile.m`
Defines the nondimensional, cross-sectional profile of the beam.  Typically taken to be a Gaussian or a tanh-tanh (smoothened top-hat) function.

## `ddn.m`
Defines the boundary value problem by a system of differential equations which takes kappa and the eigenvalue as paramters.  The system is satisfied when the parameter kappa and the associated eigenvalue are compatible.

## `odeON.m`
Ortho-normalizes the differential equation.  This step is necessary since the system is stiff and numerical error can contaminate the true solution.

## `res_eigML.m`
The objective/error function to be optimized.

## `GSBcomparison.m`
Comparison of the theoretical formulation in my PhD thesis to [Gerkema et al.](http://onlinelibrary.wiley.com/doi/10.1029/2005GL025105/full)
