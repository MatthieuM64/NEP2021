# NEP2021
Codes used in the article "Narrow escape problem in two-shell spherical domains" by M. Mangeat and H. Rieger.</br>
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
