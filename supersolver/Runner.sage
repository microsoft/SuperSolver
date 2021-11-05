from sage.all import * 
import os

#input
p=758961677
J0_base=313115042 
J0_ext=565635584
J1_base=713187954
J1_ext=136038090
supersolver=true

##########################################################
mod_poly_max = 50
##########################################################
#set up Fp and Fp2
Fp = GF(p)
c0 = -1 # tries c0 = -1,-2,... for quadratic extension
while is_square(Fp(c0)):
    c0-=1
beta = int(c0)
Fp2.<i> = GF(p^2,modulus=x^2-c0)
##########################################################

load("Solver.sage")

Solver(p, J0_base, J0_ext, J1_base, J1_ext, supersolver)
