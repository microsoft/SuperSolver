
from sage.all import *
from  random import choice
import os
import re

'''
returns the file path for the unprocessed coefficients of the l modular polynomial
'''
def FilePathModularPoly(l):
    return 'Modular_polys/phi_j_'+str(l)+'.txt'


'''
returns the file path for the preprocessed coefficients of the l modular polynomial mod p
'''
def FilePathReducedModularPoly(p,l):
    return 'reduced_modular_polys/'+str(p)+'/'+str(l)+'.txt'


'''
returns the index of the next pair of coefficient in the modular polynomial 
i.e. the indices are constructed by: looping i from 0 to l and j from 0 to i
'''
def NextIndex(i, j, l):
    if i == j:
        i += 1
        j = 0
    else:
        j +=1
    #assert (i <= l and j <= l) or (i == l+1 and (j == 0 or j == 1))
    
    return (i,j)


"""
p is a prime and L is a list of primes

for [deg1, deg2, value] in the modular polynomial file for p
write the line "value" in reduced polynoimal file
if an index does not appear in modular polynomial file, fill it with 0

this only works for prime modular polynomials because the indices only run up to (l+1,0)
"""
def ReduceModularPolynomial(p, L):

    print("start reducing modulo ", p , "of the mod polys in the followins Ls:, ", L)

    if not os.path.exists('reduced_modular_polys'):
        os.makedirs('reduced_modular_polys')
        
    if not os.path.exists('reduced_modular_polys/'+str(p)):
        os.makedirs('reduced_modular_polys/'+str(p))

    for l in L:
        # assert the modular polynomial for l exists
        assert os.path.exists(FilePathModularPoly(l)), "os.path.exists(FilePathModularPoly(l)), l="+str(l)

    for l in L:
        if os.path.exists(FilePathReducedModularPoly(p,l)):
            continue
        print("Preprocessing with ", l, p)
        with open(FilePathModularPoly(l)) as fr:
            with open(FilePathReducedModularPoly(p,l), 'w') as fw: 
                i,j = 0, 0
                for line in fr:
                    if not line:
                        continue
                    # change the string "[de1, deg2, value]" to the array [deg1, deg2, value]
                    data = [int(d) for d in re.findall(r'-?\d+', line)]
                    # write "0" if a pair of index are skipped
                    while not (i == data[0] and j== data[1]):
                        fw.write('0\n')
                        i,j = NextIndex(i,j,l)
                    fw.write(str(data[2] % p))
                    fw.write('\n')
                    i,j = NextIndex(i,j,l)
            #assert((i,j) == (l+1, 1)), "((i,j) == (l+1, 1)),"

    print("done reducing modulo p")

