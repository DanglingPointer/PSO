#### PSO
Particle Swarm Optimization algorithm implemented in C++. Input parameter values are read in from a text file. Objective and constraint functions are to be written by user.

The file `partswarmopt.h` is one of my first experiments using C++. It is very simple and doesn't always function properly. This file might be interesting for historical purposes only and should be ignored otherwise.

The second and final implementation can found in the file `pso.h` while an example of using it can be found in `example.cpp`. It allows optimization using PSO for discrete or continuous parameter ranges. In both cases objective function and constraints have to be implemented by deriving from `__interface IFitness`. Both discrete and continuous PSO can be accessed through `class Pso` which includes all necessary functions for carrying out optimization. Input format differs somewhat however:

* Discrete: all possible values for a parameter on each line.
* Continuous: minimum and maximum value for a parameter on each line.
 
In both cases the number of lines has to be equal to the number of parameters so that no newline should be present at the end of the file. Resulting output parameters will be in the same sequence as their are in the input file. There are two examples of input files named `discrete_input.txt` and `cont_input.txt` for discrete and continuous parameter ranges, respectively. 
