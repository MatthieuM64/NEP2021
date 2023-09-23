# Narrow escape problem in two-shell spherical domains
<a href="https://dx.doi.org/10.5281/zenodo.5261175"><img src="https://zenodo.org/badge/399397571.svg" alt="DOI"></a>

Codes used in the scientific publications:</br>
[1] M. Mangeat and H. Rieger, <a href="https://iopscience.iop.org/article/10.1088/1751-8121/ab4348"><i>The narrow escape problem in a circular domain with radial piecewise constant diffusivity</i></a>, J. Phys. A: Math. Theor. <b>52</b>, 424002 (2019). Preprint available on <a href='https://arxiv.org/abs/1906.06975'>arXiv</a>.</br>
[2] M. Mangeat and H. Rieger, <a href='https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.044124'><i>Narrow escape problem in two-shell spherical domains</i></a>, Phys. Rev. E <b>104</b>, 044124 (2021). Preprint available on <a href='https://arxiv.org/abs/2104.13125'>arXiv</a>.</br>
Written by M. Mangeat (2021).

## Finite Element Method (FEM)

Contains the <a href='https://freefem.org/'>FreeFem++ codes</a>. Calculate the numerical solution of the MFPT with the finite element method.</br>
Run: FreeFem++ -v 0 -nw NEP2d_twoshell_FEM.edp -parameter value.</br>
List of parameters: epsilon, Delta, D1, D2, V1, V2, Nvert, err, expT, expNE, expDELTA, expMAX (details as comments in the code).

## Kinetic Monte-Carlo (KMC)

Contains the C++ codes. Calculate the MFPT via the stochastic simulations of Brownian particles.</br>
Compile: g++ NEP2d_twoshell_KMC.cpp -lgsl -lgslcblas -lm -O3 -s -o NEP2d_twoshell_KMC.out.</br>
Run: ./NEP2d_twoshell_KMC.out -parameter=value.</br>
List of parameters: epsilon, Delta, D1, D2, V1, V2, Npart, ran (details as comments in the code).
