# To run this file:
# sage shell: load("%filename%.py")
from collections import namedtuple

load("ExtFieldSearch.sage")
load("BaseFieldSearch.sage")

'''
Edge: a tuple to store edge, corresponding field, and the prime for the isogeny
'''
Edge = namedtuple('Edge', ['J0','J1','Field','l'])

'''
checks if an Edge is valid by substituting J0, J1 into the
lth modular polynomial and asserting it be 0
'''
def CheckValidEdge(edge, p):

    X=edge.J0
    Y=edge.J1
    l = edge.l
    phi = 0

    assert os.path.exists(FilePathReducedModularPoly(p, l)), " os.path.exists(FilePathReducedModularPoly(p, l))" +str(p)+" "+str(l)
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
    assert (phi== 0), "(phi== 0)"
    return 0


'''
checks if all edges in a path is valid
'''
def CheckPathValid(listOfEdge, p):

    # check connected, assume path is nonempty
    endOfLastEdge = listOfEdge[0].J0
    for begin, end, _, _ in listOfEdge:
        assert(endOfLastEdge == begin)
        endOfLastEdge = end
    
    # check valid edge
    for edge in listOfEdge:
        CheckValidEdge(edge, p)


'''
solves the supersingular isogeny problem by 
calling Algo0Preprocess, Algo1Preprocess, Algo0, and Algo1
eventually validates the path

Also writes the solution to the problem to a file
'''
def Solver(p, J0_base, J0_ext, J1_base, J1_ext, supersolver):

    # check if directory for result is created
    if not os.path.exists('solverResults'):
        os.makedirs('solverResults')

    # check if problem solved
    filename = 'solverResults/'+str(','.join(str(num) for num in [p,J0_base,J0_ext,J1_base,J1_ext])) +'.txt'
    if os.path.exists(filename):
        print("Problem is already solved")
        return

    # bitlength of P
    size = ceil(log(p,2))
    
    # check supersingular
    assert is_prime(p), "is_prime(p)"

    _.<i>=Fp2

    initial_j0 = Fp2(J0_base + J0_ext*i)
    initial_j1 = Fp2(J1_base + J1_ext*i)

    E0 = EllipticCurve_from_j(initial_j0)
    E1 = EllipticCurve_from_j(initial_j1)

    assert E0.is_supersingular(), "E0.is_supersingular()"
    assert E1.is_supersingular(), "E1.is_supersingular()"

    #check input nodes are different
    if J0_base == J1_base and J0_ext == J1_ext:
        return [initial_j0]
    
    # algo0 preprocessing
    constants, fastest_sets, fastest_subscripts = Algo0Preprocess(p,supersolver)

    if supersolver:
        fast_ells=fastest_sets[1]
    else:
        fast_ells=[]

    print(fast_ells)

    # algo1 preprocessing
    L = Algo1Preprocess(p)

    path0, path1 = [], []

    # bothFound is true if and only if we find path from both extension
    # field element to base field element
    bothFound = True

    # algo 0 bottleneck
    if J0_ext == 0:
        path0 = [J0_base]
    else:
        path0, _, _, last_ell_1, found = Algo0(p, J0_base, J0_ext, size, constants, supersolver, fast_ells)
        bothFound = bothFound and found
        
        with open(filename, 'a') as fw:
            fw.write('Path0\n')
            fw.write(str(path0))
            fw.write('\n')

    if J1_ext == 0:
        path1 = [J1_base]
    else:
        path1, _, _,last_ell_2, found = Algo0(p, J1_base, J1_ext, size, constants, supersolver, fast_ells)
        bothFound = bothFound and found
        
        with open(filename, 'a') as fw:
            fw.write('Path1\n')
            fw.write(str(path1))
            fw.write('\n')

    if not bothFound:
        return "Unable to find"

    #algo 1 bottleneck
    base_path, l_steps = Algo1(p, path0[-1], path1[-1], L)

    # done the bottleneck computation, concatenating paths
    part1 = path0
    part2 = base_path
    part3 = path1[::-1]

    with open(filename, 'a') as fw:
        fw.write('base_path and ls\n')
        fw.write(str(base_path))
        fw.write('\n')
        fw.write(str(l_steps))
        fw.write('\n')

    path_details = []
    for index in range(len(part1)-2):
        path_details.append(Edge(part1[index],part1[index+1],'Fp2', 2))
    path_details.append(Edge(part1[len(part1)-2],part1[len(part1)-1],'Fp2', last_ell_1))
    for index in range(len(part2)-1):
        path_details.append(Edge(part2[index],part2[index+1],'Fp', l_steps[index]))
    path_details.append(Edge(part3[0],part3[1],'Fp2', last_ell_2))
    for index in range(1,len(part3)-1):
        path_details.append(Edge(part3[index],part3[index+1],'Fp2', 2))

    print("Done computing paths")
    print("------------------------")
    print("path_details:")
    
    with open(filename, 'a') as fw:  
        fw.write('path_details\n')
        for edge in path_details:
            fw.write(str(edge))
            fw.write('\n')

    for edge in path_details:
        print(edge)

    print("Start path check")
    CheckPathValid(path_details, p)
    print("Done path check")

    resulting_path = part1 + part2[1:-1] + part3
    return resulting_path
    
