# NEP2021
Codes used in the article "Narrow escape problem in two-shell spherical domains" by M. Mangeat and H. Rieger.
Written by M. Mangeat (2021).

FEM/
Contains the FreeFem++ codes (https://freefem.org/). Calculate the numerical solution of the MFPT with the finite element method.
Run: FreeFem++ -v 0 -nw $name.edp -param $param
List of parameters: epsilon, Delta, D1, D2, V1, V2, Nvert, err, expT, expNE, expDELTA, expMAX (details as comment in the code).

KMC/
Contains the C++ codes. Calculate the MFPT via the stochastic simulations of Brownian particles.
Compile: g++ name.cpp -lgsl -lgslcblas -lm -O3 -s -o name.out.
Run: ./name.out -param $param.
List of parameters: epsilon, Delta, D1, D2, V1, V2, Npart, ran (details as comment in the code).
