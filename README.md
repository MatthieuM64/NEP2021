# Narrow escape problem in two-shell spherical domains
<a href="https://dx.doi.org/10.5281/zenodo.5261175"><img src="https://zenodo.org/badge/399397571.svg" alt="DOI"></a>

Codes used in the scientific publications:</br>
[1] M. Mangeat and H. Rieger, <a href="https://doi.org/10.1088/1751-8121/ab4348"><i>The narrow escape problem in a circular domain with radial piecewise constant diffusivity</i></a>, J. Phys. A: Math. Theor. <b>52</b>, 424002 (2019).</br>
[2] M. Mangeat and H. Rieger, <a href='https://arxiv.org/abs/2104.13125'><i>Narrow escape problem in two-shell spherical domains</i></a>, submitted (2021).</br>
Written by M. Mangeat (2021).

<b>FEM/</b></br>
Contains the FreeFem++ codes (https://freefem.org/). Calculate the numerical solution of the MFPT with the finite element method.</br>
Run: FreeFem++ -v 0 -nw $name.edp -param $param</br>
List of parameters: epsilon, Delta, D1, D2, V1, V2, Nvert, err, expT, expNE, expDELTA, expMAX (details as comments in the code).

<b>KMC/</b></br>
Contains the C++ codes. Calculate the MFPT via the stochastic simulations of Brownian particles.</br>
Compile: g++ name.cpp -lgsl -lgslcblas -lm -O3 -s -o name.out.</br>
Run: ./name.out -param $param.</br>
List of parameters: epsilon, Delta, D1, D2, V1, V2, Npart, ran (details as comments in the code).
