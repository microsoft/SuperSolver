# Copyright (c) Microsoft Corporation.
# Licensed under the MIT license.

# Below functions are exactly as described in the paper

##########################################################
# Computes the degree of the level-l modular polynomial

def DegModPoly(l):

    L=list(factor(l))
    N=1
    for i in range(0,len(L)):
        N*=(L[i][0]+1)*L[i][0]^(L[i][1]-1)
    return N

##########################################################
# Computes [j,j^2,...,j^lmax] optimally 

def Poweringj(lmax,j):

    jpowers = [j,j^2] 

    if is_odd(lmax):
        for e in range(2,(lmax+1)/2+1):
            jpowers.append(jpowers[e-1]*jpowers[e-2])
            jpowers.append(jpowers[e-1]^2)
    else:
        for e in range(2,lmax/2+1):
            jpowers.append(jpowers[e-1]*jpowers[e-2])
            jpowers.append(jpowers[e-1]^2)
        jpowers.append(jpowers[lmax/2]*jpowers[lmax/2-1])

    return jpowers

##########################################################
# Evaluates \Upphi_{l,p}(X,Y) at Y=j

def EvalModPolyj(l,j,p):

    assert os.path.exists(FilePathReducedModularPoly(p, l))

    counter = [0,0,0] # records the number of base field Additions, Squares, Multiplications

    jpowers = Poweringj(DegModPoly(l),j)

    vec = [0 for i in range(0,DegModPoly(l)+1)]

    # compute phi (the modular polynomial) in X and Y
    phi = 0
    with open(FilePathReducedModularPoly(p, l)) as fr:
        a, b = 0, 0
        for line in fr:
            if line:
                data = int(line)
                #a selects power of J
                #b selects power of x 
                if a == 0: #[a,b]=[0,0] only
                    vec[0] += data;                     counter[0]+=1; 
                elif b==0: 
                    vec[0] += data * jpowers[a-1];      counter[0]+=2; counter[2]+=2 
                    vec[a] += data;                     counter[0]+=2;
                elif a == b:
                    vec[b] += jpowers[a-1] * data;      counter[0]+=2; counter[2]+=2
                else:
                    vec[b] += jpowers[a-1] * data;      counter[0]+=2; counter[2]+=2
                    vec[a] += jpowers[b-1] * data;      counter[0]+=2; counter[2]+=2
                a,b = NextIndex(a,b,l)  
    
    return vec,counter

##########################################################
# Splits into Real and Imaginary components as in the paper

def RealAndImaginary(vec):

    real=[]
    imag=[]

    for v in vec:
        vsplit=Fp2(v).polynomial().list()
        if len(vsplit) == 0:
            real.append(Fp2(0))
            imag.append(Fp2(0))
        elif len(vsplit)==1:
            real.append(vsplit[0])
            imag.append(Fp2(0))
        else:
            real.append(vsplit[0])
            imag.append(vsplit[1])    
    
    assert imag[len(imag)-1] == 0
    imag.pop(len(imag)-1)

    return real,imag

##########################################################
# Computes the inversion-free polynomial GCD as in the paper

def InvFreeGCD(g,h):

    counter = [0,0,0] # records the number of base field Additions, Squares, Multiplications

    r,s=g,h

    lcr,lcs=r[len(r)-1],s[len(s)-1]
    r,s=[lcs*x for x in r],[lcr*x for x in s];          counter[2]+=len(r)+len(s)

    while len(r) > 2 and r != s:
        degDiff=len(r)-len(s)
        for i in range(len(s)):                         
            r[i+degDiff]=r[i+degDiff]-s[i];             counter[0]+=1
        while r[len(r)-1] == 0:
            r.pop(len(r)-1)
        lcr,lcs=r[len(r)-1],s[len(s)-1]
        r,s=[lcs*x for x in r],[lcr*x for x in s];      counter[2]+=len(r)+len(s)
        if len(r) < len(s):
            r,s=s,r
    
    if len(r) == 2 and r != s:
        return false, counter
    else:
        return true, counter

##########################################################

def NeighbourInFp(l,j,p):

    vec,count1=EvalModPolyj(l,j,p)
    real,imag=RealAndImaginary(vec)
    bool,count2= InvFreeGCD(real,imag)
    count=[count1[0]+count2[0],count1[1]+count2[1],count1[2]+count2[2]]

    return bool,count,vec

##########################################################

