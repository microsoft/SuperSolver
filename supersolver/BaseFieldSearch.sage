# To run this file:
# sage shell: load("%filename%.sage")

'''
Given a prime p and two field elements J0, J1 in Fp
find a path between them (Delfs-Galbraith Algorithm)
'''

load("ModularPolynomials.sage")

"""
given p, find the primes less than B that has kronecker symbol 1
then, write the coefficients of lth modular polynomial
mod p to preprocessed files
"""
def Algo1Preprocess(p):
    
    assert is_prime(p)
    d = -4 * p
    # B = math.floor(6*(math.log(abs(d))**2)) #Delfs-Galbraith upper bound
    B = mod_poly_max 

    L = [l  for l in prime_range(1, B) if kronecker_symbol(-p, l)==1 and l <= B]
    return L

'''
Delfs-Galbraith Algorithm: find a path between J0 and J1
'''
def Algo1(p, J0, J1, L):

    x = SupersingularModule(p)

    X = PolynomialRing(Fp, 'x').gen()

    l_step = [[],[]]
    l_step_return = []
    S = [[Fp(J0)], [Fp(J1)]]

    # O(1) lookup to see if elemnt is in S
    Slookup = [set([Fp(J0)]), set([Fp(J1)])] 

    disjoint = True
    i = 0

    while disjoint:
        l = 0
        roots = []
        while len(roots) == 0:
            l = choice(L)
            Y = S[i][-1]

            # assumes ReduceModularPolynomial(p,L) is already called and l is in L
            assert os.path.exists(FilePathReducedModularPoly(p, l))

            # compute phi (the modular polynomial) in X and Y
            phi = 0
            with open(FilePathReducedModularPoly(p, l)) as fr:
                a, b = 0, 0
                for line in fr:
                    if line:
                        data = int(line)
                        if a == b:
                            phi += X^a * Y^b * data
                        else:
                            phi += X^a* Y^b * data + X^b * Y^a * data
                        a,b = NextIndex(a,b,l)
            roots = [root[0] for root in phi.roots() if root[0] not in Slookup[i]]
            # now, if all the roots in roots is seen previously on the same side then roots would be
            # empty and we pick a new l

        # j is assigned a random root in the modular polynomial
        j = choice(roots)
        l_step[i].append(l)
        S[i].append(j)
        Slookup[i].add(j)

        if j in Slookup[1-i]:
            #print("found!", j)
            disjoint = False
            index0 = S[0].index(j)
            index1 = S[1].index(j)

            # finds the position of the common vertex, and concatenate corresponding 
            # part of the paths
            S = S[0][:index0+1] + (S[1][:index1])[::-1]
            l_step_return = l_step[0][:index0] + (l_step[1][:index1])[::-1]
        
        i = 1 - i
    
    print('------------------')

    return [S, l_step_return]
