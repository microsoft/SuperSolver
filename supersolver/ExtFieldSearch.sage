# Copyright (c) Microsoft Corporation.
# Licensed under the MIT license.

'''
Given a prime p and a field element (base + ext*i) in Fp^2
find a path from the field element to a base field element.
'''

load('NeighbourInFp.sage')
load('ModularPolynomials.sage')

'''
Returns the appropriate tree depth for 
the iterateive graph search in Algo0
'''
def RecursionDepth(p): 
    #return ceil(log(p,2)) # gives precision errors
    return len(bin(p))-2   # default depth is the approximate diameter of the graph

'''
Does a plain square and multiply. Users obtaining precise data for a fixed 
prime p can modify this function using tailored exponentiation methods
which should slightly improve the cost of the 2-isogeny step.
'''
def Expon(a,n):

    counter=[0,0,0]
    nbits = [int(x) for x in bin(n)[2:]]    
    x=a
    for i in range(1,len(nbits)):  
        x=x^2;                          counter[1] += 1
        if nbits[i] == 1:
            x*=a;                       counter[2] += 1
    return [x,counter]

'''
The functions below for computing optimised square roots in GF(p^2)
follow the methods described in Michael Scott's "Tricks of the trade" 
paper. They implement the "progenitor" method with the Tonelli-Shanks algorithm
'''

def TonelliShanksUpdate(x,y,consts):

    counter=[0,0,0]
    a0= int(consts[0])
    a1= consts[1]
    a2= consts[2]
    a3= consts[3]
    s0=y
    s1=x
    for _ in range(2, a0+1):
        s0=s0^2;                        counter[1]+=1
        s1=s1^2;                        counter[1]+=1
    s0=s0^2;                            counter[1]+=1
    s0*=s1;                             counter[2]+=1
    if s0 != 1:
        x*=a2;                          counter[2]+=1
        y*=a3;                          counter[2]+=1
    s1=x*y;                             counter[2]+=1
    s2=s1*y;                            counter[2]+=1
    for k in range(a0, 1, -1):
        s3=s2
        for i in range(3,k+1):
            s3*=s3;                     counter[2]+=1
        if s3 != 1:
            s1*=a1;                     counter[2]+=1
        a1*=a1;                         counter[2]+=1
        if s3 != 1:
            s2*=a1;                     counter[2]+=1
    return [s1,s0,counter]

# computes the square root of an extension field element
# increments the operations counter accordingly

def Fp2Sqrt(element,constants):

    c0 = constants[0]
    c1 = constants[1]
    c2 = constants[2]
    c3 = constants[3]
    c4 = constants[4]
    c5 = constants[5]
    c6 = constants[6]
    c7 = constants[7]
    c8 = constants[8]

    counter=[0,0,0]
    tuple = element.polynomial().list()
    t0=tuple[0]
    t1=tuple[1]

    t2=t0^2;                            counter[1]+=1
    t3=t1^2;                            counter[1]+=1
    t3=c0*t3;                           counter[2]+=1
    t2=t2-t3;                           counter[0]+=1

    t3,exp_count=Expon(t2,c4)
    counter[0]+=exp_count[0]; counter[1]+=exp_count[1]; counter[2]+=exp_count[2]

    t2,t3,upd_count=TonelliShanksUpdate(t2,t3,[c2,c3,c0,c5])
    counter[0]+=upd_count[0]; counter[1]+=upd_count[1]; counter[2]+=upd_count[2]

    t3=t0+t2 ;                          counter[0]+=1
    t3=t3*c1 ;                          counter[2]+=1
    t0,exp_count=Expon(t3,c4)
    counter[0]+=exp_count[0]; counter[1]+=exp_count[1]; counter[2]+=exp_count[2]
    t2=t0

    for _ in range(1,c7+1):
        t2=t2^2  ;                      counter[1]+=1
    
    t2*=c1;                             counter[2]+=1
    t2*=t1;                             counter[2]+=1
    t1,exp_count= Expon(t3,c6)
    counter[0]+=exp_count[0]; counter[1]+=exp_count[1]; counter[2]+=exp_count[2]
    t2*=t1;                             counter[2]+=1
    t3,t0,upd_count=TonelliShanksUpdate(t3,t0,[c2,c3,c0,c5])
    counter[0]+=upd_count[0]; counter[1]+=upd_count[1]; counter[2]+=upd_count[2]
    t2*=t3;                             counter[2]+=1

    if t0 == 1:
        t0=t3
        t3=t2
    else:
        t0=t2
        t3*=c8;                         counter[2]+=1

    return [t0,t3,counter]


# given J in 2-isogeny graph and J0 one of its neighbours
# generate the two other neighbours using the quadratic
# function derived from modular polynomials
# increments the operations counter accordingly

def GenerateValidNeighbourQuadratic(J0, J,constants):

    j2=J-1488      #A
    j2*=J          #3M+5A
    u0=j2+162000   #A
    u1=J0+J0       #2A   
    u1+=j2         #2A
    u1-=54000      #A
    u1*=u0         #3M+5A
    j2*=70016      #2M
    j1=J*57316692  #2M
    j2-=j1         #2A
    j1=J0^2        #3M+5A
    j2-=j1         #2A
    j1=j2+j2       #2A
    u1+=j1         #2A
    u1+=j2         #2A
    j2=u0-J0       #2A

    _.<i>=Fp2
    counter=[0,0,0]
    
    if u1 not in base_field(Fp2):
        t0,t3,counter=Fp2Sqrt(u1,constants)
        u1=t0+t3*i
    else:
        u1=base_field(Fp2)(u1).sqrt()
    
    j1=j2+u1                 #2A
    j1*=constants[1]         #3M+5A
    j2-=j1                   #2A 
                             #16M+43A
    counter[2] += 16
    counter[0] += 43

    return [j1,j2,counter]


'''
Finds path in the 2-isogeny graph from base + ext*i
to a vertex in base field using a binary tree and a running
stack. See the paper for further details. 
'''

def Algo0(p, base, ext, size, constants, supersolver, fast_ells):

    _.<i>=Fp2

    beta = constants[-1]
    c0 = constants[0]
    c1 = constants[1]
    c2 = constants[2]
    c3 = constants[3]
    c4 = constants[4]
    c5 = constants[5]
    c6 = constants[6]
    c7 = constants[7]
    c8 = constants[8]

    # bitlength of p
    #size = ceil(log(p,2)) #gives precision errors
    size = len(bin(p))-2

    initial_j = base + ext*i

    E = EllipticCurve_from_j(initial_j)
    assert not initial_j in GF(p), "initial_j not in GF(p)"

    counter = [0,0,0] # records the number of base field Additions, Squares, Multiplications
    nodesVisited = 1 # records total number of nodes visited

    # use modular polynomials to generate the three
    # neighbours of currentJ. This function is only called
    # once in each function call

    def generateValidNeighbourModuli(currentJ):
        listOfNext = []
        E=EllipticCurve_from_j(currentJ)
        for P in E(0).division_points(2):
            phi = E.isogeny(P)
            j = phi.codomain().j_invariant()
            if j!= currentJ:
                listOfNext.append(j)
        return listOfNext


    # the following algorithm is an iterative DFS to
    # search a complete binary tree of depth size
    # starting at initial_j

    stack=[]
    path=[initial_j] 
    right_branch=[]
    stack.append(0)
    found= False
    index = 1

    #to store the isogeny degree that finds the subfield node
    last_ell = 2

    res = generateValidNeighbourModuli(initial_j)
    if res[1] in Fp:
        path.append(res[1])
    else:
        path.append(res[0])
    right_branch.append(res[1])
    index += 1
    stack.append(0)

    if path[len(path)-1] in Fp:
        found= True

    percentage = 0
    sqrt_p = int(math.sqrt(p))

    # boolean flag to print progress
    print_progress = RecursionDepth(p) >= 40
    while len(path) != 1 and not found:

        if index == size + 2:

            path.pop()
            index -= 1
            stack[-1] += 1 
            right_branch.pop()

        elif index == len(stack) :

            if supersolver:
                j = path[len(path)-1]
                if fast_ells != []:
                    for ell in fast_ells:
                        if not found:
                            found,count,vec=NeighbourInFp(ell,j,p)
                            counter[0]+=count[0]; counter[1]+=count[1]; counter[2]+=count[2];
                            if found:
                                last_ell=ell
                                X = PolynomialRing(Fp2, 'X').gen()
                                f = 0
                                for i in range(len(vec)):
                                    f+=vec[i]*X^i
                                rts=f.roots()
                                for i in rts:
                                    if i[1] in Fp:
                                        path.append(i[0])
                                        break
                                break

            if not found:
                res = GenerateValidNeighbourQuadratic(path[len(path)-2],path[len(path)-1],constants)
                counter[0]+=res[2][0]; counter[1]+=res[2][1]; counter[2]+=res[2][2]; 
                nodesVisited += 2
                if res[1] in Fp:
                    path.append(res[1])
                    print(res[0])
                else:
                    path.append(res[0])
                right_branch.append(res[1])
                index += 1
                stack.append(0)

                if path[len(path)-1] in Fp:
                    found= True

        elif stack[index] == 0:
            path.append(right_branch.pop())
            stack[index]= 1
            index += 1
            if path[len(path)-1] in Fp:
                found= True
        else: 
            path.pop()
            stack.pop()
            index -= 1

        if print_progress:
            new_percentage = (100 * nodesVisited) // sqrt_p
            if new_percentage > percentage:
                percentage = new_percentage
                print("% done (expected):", percentage, ", nodes visited/expected:", nodesVisited,"/",sqrt_p )

    if supersolver:
        if fast_ells != []:
            counter[2]+=fast_ells[-1] #for the j-powering

    return path, counter, nodesVisited, last_ell, found

'''
the preprocessing steps for Algo0,
compute some Tonelli-Shanks constants needed for Algo0
'''
def Algo0Preprocess(p,supersolver):

    ReduceModularPolynomial(p, [l for l in range(2,mod_poly_max+1)]) 

    C0=Fp(c0) # smallest QNR in GF(p)
    c1=Fp(1/2) # constant for division by 2  
    c2=1 # Tonelli-Shanks "e" in Michael Scott's paper
    while (p-1) % 2^(c2+1) == 0:
        c2+=1
    c3=C0^((p-1)/(2^c2)) # "z" in Michael Scott's paper
    c4=(p-1-2^c2)/(2^(c2+1)) # y=x^c4 in Michael Scott's paper
    c4=Fp(c4)
    c6=2^c2-1
    c7=c2+1
    c8=1/C0 #1/c0 in Michael Scott's paper
    c5=C0**int(c4)
    constants=[C0,c1,c2,c3,c4,c5,c6,c7,c8,beta]
    fast_ells=[]
    ell_counts=[]
    fastest_sets=[[]]
    fastest_subscripts=[]

    def sort_list(list1, list2):
 
        zipped_pairs = zip(list2, list1)
        z = [x for _, x in sorted(zipped_pairs)]
        return z

    if supersolver:

        j1,j2,counter=GenerateValidNeighbourQuadratic(Fp2.random_element(),Fp2.random_element(),constants)
        muls_step=int((counter[1]+counter[2])/2) 

        print(" ")
        print("One step in the 2-isogeny graph takes ~ ", muls_step, "Fp-multiplications ")
        print(" ")
        ell=1
        
        while ell < mod_poly_max-1:
            ell=ell+2
            _,count,_=NeighbourInFp(ell,Fp2.random_element(),p)
            muls_ell=int(count[1]+count[2])
            if muls_ell < ((ell+1)*muls_step):
                fast_ells.append(ell)
                ell_counts.append(muls_ell)
                print("Rapid inspection of ", DegModPoly(ell), " nodes in the ", ell, "-isogeny graph takes at most ", round(muls_ell/(ell+1),2), "Fp-multiplications per node")   
    
        bits_array=[]
        av_cost=[]

        # exhaustively optimise the cost of inspecting sets of \ell_i

        for i in range(2^len(fast_ells)):
            bits=bin(i)
            bits=bits[2:]
            while len(bits) < len(fast_ells):
                bits='0'+bits
            bits_array.append(bits)
            sum_numer=2*muls_step
            sum_denom=2
            for j in range(len(bits)):
                if bits[j] == '1':
                    sum_numer+=ell_counts[j]
                    sum_denom+=fast_ells[j]+1
            av_cost.append(round(sum_numer/sum_denom,4))

        av_cost_sorted=sorted(av_cost)
        bits_array_sorted=sort_list(bits_array,av_cost)

        #constructing new sets of fast_ells based on the minimum cost

        for j in range(5): #number of fast_sets to look at
            t = []
            for i in range(len(fast_ells)):
                if bits_array_sorted[j][i] == '1':
                    t.append(fast_ells[i])
            fastest_sets.append(t)
            fastest_subscripts.append(int(bits_array_sorted[j][::-1],2))

        print(" ")
        
    return constants,fastest_sets,fastest_subscripts

