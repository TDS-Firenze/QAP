# The Quadratic Assignment Problem (QAP)

In this repository there are codes connected to my Master's Thesis (https://github.com/Tommaso-Mannelli-Mazzoli/masters-thesis).
I studied the Quadratic Assignment Problem (QAP) and I implemented heuristic and metaheuristic algorithm to find numerical "good" solutions.

I used Fortran 90 with compiler  Silverfrost ftn95.


The files are:

* QAP: The main program. 
* Greedy1: this is an implementation of Greedy1 algorithm;
* Greedy2: this is an implementation of Greedy2 algorithm;
* Greedy3: this is an implementation of Greedy3 algorithm;
* Local_Search: contains the definition of Δ, Δ¹, Δ²  functions and the following algorithms:
  - 2-optimum: first improvement;
  - 2-optimum: best improvement;
  - 3-optimum: first improvement;
  - 3-optimum: best improvement.
* ACO: this is an implementation of Ant Colony Optimization algorithm;
* VNS: this is an implementation of GVNSfirst and GVNSbest algorithms;
* TABU_SEARCH: this is an implementation of Tabu Search Algorithm
* todos: text file that contains the name of the instances used to read the matrices A and B;
* Todos_SLN: text file that contains the name of the instances used to read the Best Known Solution (BKS).
