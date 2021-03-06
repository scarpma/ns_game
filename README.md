# Incompressible Navier-Stokes equation solver

Fortran implemetation of Incompressible Navier-Stokes equations using *Jos Stam Algorithm*. This program was developed as project for a master degree (theoretical Physics) exam "Fisica Computazionale".

This algorithm is a popular one in computer graphics, but does not have good numeric properties (convergence and stability).

The aim of this program is to implement and test these limitations.

Currently, the repository is splitted into two branches, the general *real-space* version, and the *Fourier-space* one, restricted to *periodic boundary conditions*. In the real-space version it is easier to implement generic boundary conditions. This is proved through the implementation of a *round obstacle inside a wind-tunnel*, as well as many others. For the Fourier space, on the contrary, only periodic boundary conditions are possible (since the *FFT algorithm is used*, i.e. the fast discrete Fourier transform). For this reason we test the case with a stocastic forcing mechanism in order (to try) to reach a stable turbulence regime in which there is a full energy transfer towards small scales (further analysis code is needed to analize this).

Numerical checks are made (for the real-space version) with respect to high-precision numerical solutions provided by Ghia et al. in their papers.
