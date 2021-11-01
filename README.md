# Wall-Models-for-LES
This repository contains MATLAB implementations of the following wall models (WM) for incompressible LES solvers. 

1. ODE equilibrium WM (traditional finite volume implementation)
2. ODE equilibrium WM (new spectral implementation)
3. Integral WM

For ODE equilibrium WM, the spectral implementation should be a preferred choice. Compared to the finite-volume implementation, it reduces the wall-modeling cost significantly along with dramatic improvement in parallel efficiency when run on multiple cores. 

Channel flow DNS data from Johns Hopkins Turbulence Database (JHTDB) and UT-Austin database are included for offline a priori validation of these wall models.

NOTE: ArXiv link to the paper describing these wall models will be posted below soon. Please check back later.
