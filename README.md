# SuperSolver

SuperSolver is a Sage/Python library that solves the general supersingular isogeny problem using the Delfs-Galbraith algorithm. It accompanies the paper “SuperSolver: accelerating the Delfs-Galbraith algorithm with fast subfield root detection”, by Maria Corte-Real Santos, Craig Costello and Jia Shi. 

# Running the code

When all of the files are in the Sage directory, the command 

load("Runner.sage")

will perform the Delfs-Galbraith algorithm to find the isogeny between E1/GF(p^2) and E2/GF(p^2). There are default parameters in Runner.sage that can be modified to any instance of the supersingular isogeny problem. The 6 inputs specified by the user are 

p: the characteristic of the field GF(p^2)
j10: an element of GF(p)
j11: an element of GF(p)
j20: an element of GF(p)
j21: an element of GF(p)
supersolver: a boolean

The j-invariant of E1/GF(p^2) is j1=j10+j11.alpha and the j-invariant of E2/GF(p^2) is j2=j20+j21.alpha, where GF(p^2)=Fp(alpha) is the quadratic extension field specified in Runner.sage. Tiny primes (i.e. p<100) should be avoided since there might not be any supersingular j-invariants in GF(p^2)\GF(p), and the code assumes both input j-invariants are in GF(p^2); note that the supersingular isogeny problem is trivial for such small primes. 

The supersolver boolean variable is true by default. Disabling it will call the original Delfs-Galbraith walk in the 2-isogeny graph, without the SuperSolver optimisations described in the paper above. 

# Bugs, Questions, Comments? 

Please contact Maria Corte-Real Santos (maria.santos.20@ucl.ac.uk) and/or Craig Costello (craigco@microsoft.com) and/or Jia Shi (janeshi99@gmail.com). 
