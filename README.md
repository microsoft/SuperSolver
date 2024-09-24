# SuperSolver

SuperSolver is an algorithm written in Sage/Python that solves the general supersingular isogeny problem using the Delfs-Galbraith algorithm. It accompanies the paper "[Accelerating the Delfs–Galbraith Algorithm with Fast Subfield Root Detection](https://link.springer.com/chapter/10.1007/978-3-031-15982-4_10)”, by Maria Corte-Real Santos, Craig Costello and Jia Shi. 

# Downloading the Modular Polynomial files

The files containing the modular polynomials (in [supersolver/Modular_polys/](https://github.com/microsoft/SuperSolver/tree/main/supersolver/Modular_polys)) are stored using Git Large File Storage (LFS). For the code to run properly, first install `git-lfs`. Then, after navigating to the folder containing this repository, run `git lfs pull`.

The modular polynomials can also be accessed from [Sutherland's database](https://math.mit.edu/~drew/ClassicalModPolys.html). However, we remark that we fix our ordering so the first line always contains the coefficient corresponding to `[0,0]`. This means that for some $N$, the modular polynomial of level $N$ is stored in reverse order compared to Sutherland's database. For example, Sutherland stores the modular polynomial $\Phi_2(X,Y)$ of level $2$ as:

```
[3,0] 1
[2,0] -162000
[2,1] 1488
[2,2] -1
[1,0] 8748000000
[1,1] 40773375
[0,0] -157464000000000
```
We store this polynomial as:

```
[0,0] -157464000000000
[1,0] 8748000000
[1,1] 40773375
[2,0] -162000
[2,1] 1488
[2,2] -1
[3,0] 1
```

# Running the code

When all of the files are in the Sage directory, the command 

load("Runner.sage")

will perform the Delfs-Galbraith algorithm to find the isogeny between E1/GF(p^2) and E2/GF(p^2). There are default parameters in Runner.sage that can be modified to any instance of the supersingular isogeny problem. The code will generate new folders in the local directory for the reduced modular polynomials and solution. 

The 6 inputs specified by the user are 

p: the characteristic of the field GF(p^2)

j10: an element of GF(p)

j11: an element of GF(p)

j20: an element of GF(p)

j21: an element of GF(p)

supersolver: a boolean

The j-invariant of E1/GF(p^2) is j1=j10+j11.alpha and the j-invariant of E2/GF(p^2) is j2=j20+j21.alpha, where GF(p^2)=Fp(alpha) is the quadratic extension field specified in Runner.sage. Tiny primes (i.e. p<100) should be avoided since there might not be any supersingular j-invariants in GF(p^2)\GF(p), and the code assumes both input j-invariants are in GF(p^2); note that the supersingular isogeny problem is trivial for such small primes. 

The supersolver boolean variable is true by default. Disabling it will call the original Delfs-Galbraith walk in the 2-isogeny graph, without the SuperSolver optimisations described in the paper above. 

The file ExperimentAlgo0.sage can be run to obtain data on many instances of the problem for a given prime. After specifying a prime in ExperimentAlgo0.sage, the command

load("ExperimentAlgo0.sage") 

followed by the command 

Experiment(p,N) 

will run N (an integer) instances of the subfield search problem for pseudo-random supersingular j-invariants over GF(p^2), as in the paper.

# Bugs, Questions, Comments? 

Please contact Maria Corte-Real Santos (maria.santos.20@ucl.ac.uk) and/or Craig Costello (craigco@microsoft.com) and/or Jia Shi (janeshi99@gmail.com). 

# Code of Conduct

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/). For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or contact opencode@microsoft.com with any additional questions or comments.

## Contributing

This project welcomes contributions and suggestions.  Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit https://cla.opensource.microsoft.com.

When you submit a pull request, a CLA bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., status check, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.

# License

Copyright (c) Microsoft Corporation. All rights reserved.

Licensed under the [MIT](https://github.com/microsoft/vscode/blob/main/LICENSE.txt) license.
