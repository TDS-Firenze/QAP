# The Quadratic Assignment Problem (QAP)

In this repository there are codes connected to my Master's [Thesis](https://github.com/Tommaso-Mannelli-Mazzoli/masters-thesis).

I studied the Quadratic Assignment Problem (QAP) and I implemented heuristic and metaheuristic algorithm to find numerical "good" solutions.

I used Fortran 90 with compiler  Silverfrost ftn95.


The files are:

* **QAP**: A complete file that contains every function and subroutine;
* **main**: The main program.
## Subroutines

* **Greedy1**: this is an implementation of Greedy1 algorithm;
* **Greedy2**: this is an implementation of Greedy2 algorithm;
* **Greedy3**: this is an implementation of Greedy3 algorithm;
* **Local_Search**: contains the following algorithms:
  - 2-optimum: first improvement;
  - 2-optimum: best improvement;
  - 3-optimum: first improvement;
  - 3-optimum: best improvement.
* **ACO**: this is an implementation of Ant Colony Optimization algorithm;
* **VNS**: this is an implementation of GVNSfirst and GVNSbest algorithms;
* **TABU**: this is an implementation of Tabu Search Algorithm.

## Functions
* **funobj**: evaluates the objective function value;
* **delta**:  evaluates the difference of Objective Function Value after the 2-exchange {i1,i2} --> {i2,i1};
* **delta1**: evaluates the difference of Objective Function Value after the 3-exchange {i1,i2,i3} --> {i2,i3,i1};
* **delta2**: evaluates the difference of Objective Function Value after the 3-exchange {i1,i2,i3} --> {i3,i1,i2}.

## Instancies
This folder contains the Problem Instances and Solutions of QAP. The files here are of two types:
* **.dat**: the data of the instance.

 The format of the problem data is:
 
>*s*     *n* 


>***A***


>***B***

where *s* = 1 if the problem is symmetric (0 if it is not), *n* is the size of the instance, and ***A*** and ***B*** are flow and distance matrix.
* **.sln**: the Best Known Solution of the instance.
The format of these files is

>*n*    *sol* 

>*p*


where *n* gives the size of the instance, *sol* is the objective function value and *p* a corresponding permutation.

## Other files
* **todos**: text file that contains the name of the instances used to read the matrices A and B.
* **Todos_SLN**: text file that contains the name of the instances used to read the Best Known Solution (BKS).
