from sage.all import * 
import os

# Input parameters that define the supersingular isogeny problem

p=758961677
j10 =313115042 
j11=565635584
j20=713187954
j21=136038090
supersolver=true

##########################################################
# Modular polynomials currently restricted to level 50
# This can be increased for larger instances, but the modular
# polynomials (taken from Sutherland's database) 
# must then be inserted into the Modular_polys folder

mod_poly_max = 50

##########################################################
# set up Fp and Fp2
Fp = GF(p)
c0 = -1 # tries c0 = -1,-2,... for quadratic extension
while is_square(Fp(c0)):
    c0-=1
beta = int(c0)
Fp2.<i> = GF(p^2,modulus=x^2-c0)
##########################################################

load("Solver.sage")

Solver(p, j10, j11, j20, j21, supersolver)
