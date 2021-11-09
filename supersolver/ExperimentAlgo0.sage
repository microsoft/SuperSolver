# Copyright (c) Microsoft Corporation.
# Licensed under the MIT license.

load("ExtFieldSearch.sage")
import os

'''
Given bitsize, numberOfPrimes, numberOfinstances,
generate $numbersOfPrimes$ primes of appropriate bitlength
and for each prime, generate $numberOfinstances$ random 
extension field elements and perform Algo0 on each instance.
returns statistics
'''

p=previous_prime(2^28) #can be changed by user for experiments

Fp = GF(p)
c0 = -1 # tries c0 = -1,-2,... for quadratic extension
while is_square(Fp(c0)):
    c0-=1
Fp2.<i> = GF(p^2,modulus=x^2-c0)

'''
given p and beta (which is the modulus), generate a random j-invariant
'''
def Getj(p):

    j_inv=supersingular_j(Fp2)
    E=EllipticCurve_from_j(j_inv)

    if E.cardinality() % 3 != 0:
        E = E.quadratic_twist()
        assert E.cardinality() % 3 == 0 and E.is_supersingular()

    for _ in range(20):
        P = choice(E(0).division_points(3))
        phi = E.isogeny(P)
        E = phi.codomain()
    
    while E.j_invariant() in Fp:
        for _ in range(20):
            P = choice(E(0).division_points(3))
            phi = E.isogeny(P)
            E = phi.codomain()
    assert E.is_supersingular()

    result = [Fp(num) for num in E.j_invariant().polynomial().list()]
    return result


'''
generate $numberOfPrimes$ number of primes of bitsize
'''
def GeneratePrimesGivenBitsize(bitsize, numberOfPrimes):
    assert bitsize > 2, "bitsize > 2"
    # returns a size numberOfPrimes list of primes between 2^bitsize-1, 2^(bitsize-1)
    primes = set()

    # remove the below line for primes with bitlength larger than 60, as it could take very long time
    assert prime_pi(2^bitsize - 1) - prime_pi(2^(bitsize-1)) >= numberOfPrimes, "not enough primes in given range"

    while len(primes) < numberOfPrimes:
        p = random_prime(2^bitsize-1, false, 2^(bitsize - 1))
        assert is_prime(p), "is_prime(p)"
        primes.add(p)

    assert len(primes) == numberOfPrimes

    return primes


'''
generates $numberOfinstances$ j-invariants for prime
'''
def GeneratejsGivenPrime(prime, numberOfinstances, beta):
    j_invs = set()
    while len(j_invs) < numberOfinstances:
        base, ext = Getj(prime)
        if ext != 0:
            j_invs.add((base, ext))

    assert len(j_invs) == numberOfinstances, "len(j_inv) == numberOfinstances"

    return list(j_invs)

'''
Creates the directory and the files to write solutions and stats
If file with the same bitlength already exists, erase it.
'''
def PrepareExperimentFiles(bitlength,index):

    if not os.path.exists('ExperimentResults'):
        os.makedirs('ExperimentResults')

    if not os.path.exists('ExperimentResults/'+str(bitlength)+'-'+str(index)):
        os.makedirs('ExperimentResults/'+str(bitlength)+'-'+str(index))
    
    filepathSolutions = 'ExperimentResults/'+str(bitlength)+'-'+str(index)+'/solutions.txt'
    filepathStats = 'ExperimentResults/'+str(bitlength)+'-'+str(index)+'/stats.txt'

    with open(filepathSolutions, 'w') as fw:
        fw.truncate()
    with open(filepathStats, 'w') as fw:
        fw.truncate()


def WriteSolution(bitlength, p, j_inv, index, fast_ells, path):
    filepath = 'ExperimentResults/'+str(bitlength)+'-'+str(index)+'/solutions.txt'
    assert os.path.exists(filepath)
    with open(filepath, 'a') as fw:
        fw.write('p='+str(p)+', j_inv='+str(j_inv)+', path='+str(path)+', fast ells index='+str(index)+' ,fast_ells='+str(fast_ells)+ '\n')
    return

def WriteStats(bitlength,index, totalInstances, avgM, avgS, avgA, avgNodes, maxNodes, final=False): 
    filepath = 'ExperimentResults/'+str(bitlength)+'-'+str(index)+'/stats.txt'
    assert os.path.exists(filepath)
    with open(filepath, 'a') as fw:
        if final:
            fw.write('FINAL RESULT\n')
        fw.write('totalInstances='+str(totalInstances)+'\n')
        fw.write('avgM='+str(int(avgM+avgS))+'\n')
        fw.write('avgNodesVisited='+str(avgNodes)+'\n')
    return

'''
the full Experiment
'''

def Experiment(p, numberOfinstances):

    supersolver=1
    bitsize=ceil(log(p)/log(2))

    constants,fastest_sets, fastest_subscripts = Algo0Preprocess(p,supersolver)
    instances = GeneratejsGivenPrime(p, numberOfinstances, constants[-1])

    depth = RecursionDepth(p)

    for fast_ells in fastest_sets:

        totalCounter, totalNodesVisited = [0,0,0], 0
        totalCounterSquared, totalNodesVisitedSquared = [0,0,0], 0
        totalInstances = 0
        minNodesVisited, maxNodesVisited = float("inf"), float("-inf")
        PrepareExperimentFiles(bitsize,fastest_sets.index(fast_ells))

        for j_inv in instances:

            # if unable to find basefield node, reduce number of instances
            # and do not update the statistics with this instance

            path, counter, nodesVisited, last_ell, found = Algo0(p, j_inv[0], j_inv[1], depth, constants, supersolver,fast_ells)
                                                        
            if found:   
                WriteSolution(bitsize, p, j_inv,fastest_sets.index(fast_ells),fast_ells, path)
                totalInstances += 1
                totalCounter[0] += counter[0]
                totalCounter[1] += counter[1]
                totalCounter[2] += counter[2]
                totalNodesVisited += nodesVisited

                if nodesVisited < minNodesVisited:
                    minNodesVisited = nodesVisited
                if nodesVisited > maxNodesVisited:
                    maxNodesVisited = nodesVisited

        averageCounter = [round(i / totalInstances,2) for i in totalCounter]
        averageNodesVisited = int(totalNodesVisited / totalInstances)

        WriteStats(bitsize, fastest_sets.index(fast_ells), totalInstances, averageCounter[2], averageCounter[1], averageCounter[0], averageNodesVisited, maxNodesVisited, True) 

        print("Bitsize:",bitsize, " with ",numberOfinstances," instances:" )

        resultMean = [round(frac,2) for frac in averageCounter + [averageNodesVisited]]

        print(totalInstances, " number of instances")
        print("Mean M:", resultMean[2])
        print("Mean S:", resultMean[1])
        print("Mean A:", resultMean[0])

        print("Average nodes visited:", resultMean[3])
        print("Max nodes visited:", maxNodesVisited)
        print("Min nodes visited:", minNodesVisited)
        print("----------------------------------------")

    return "done"
