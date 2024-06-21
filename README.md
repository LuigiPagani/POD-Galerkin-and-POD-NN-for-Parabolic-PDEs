## Overview
This repository contains the implementation of a 2D unsteady convection-diffusion simulation, primarily focused on modeling the dispersion of pollutants in water influenced by turbines. The simulation employs FEniCS for solving PDEs and implements a reduced basis method through the POD-Galerkin scheme for efficient computation.

## Features
- **Convection-Diffusion Solver:** Simulates the spread of a chemical pollutant using the fluid velocities obtained from the Stokes' solver.
- **Reduced Basis Method:** Implements the POD-Galerkin scheme to reduce computational overhead while maintaining accuracy, and performing the time stepping in the reduced basis space.
- **DL-ROM** Implements the DL_ROM scheme to reduce computational overhead while maintaining accuracy, and treating time as a parameter, allowing to query the solver indipendently at each time step.
- **Performance Analysis:** Evaluates the efficiency and accuracy of the reduced model compared to the full-order model (FOM).


### Configuration
Adjust the simulation parameters, including the strength coefficients of the turbines and the final simulation time.

