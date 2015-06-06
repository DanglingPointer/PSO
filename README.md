#### PSO
Particle Swarm Optimization algorithm implemented in C++. Input parameter values are read in from a text file. Objective and constraint functions are to be written by user as part of a class derived from `IFitness`.

The file `partswarmopt.h` is one of my first experiments using C++. It is very simple and doesn't always function properly. The function calculating Inertia weight is a stub. This file is only interesting for historical purposes and should be ignored otherwise.

The second and final implementation is more advanced and performs well. It can found in the file `pso.h` while an example of using it can be found in `Example_discrete.cpp`. User has to provide a text file with discret input parameter values. All possible values for a parameter has to be on the same line. Hence, the number of lines should be equal to the number of parameters. Objective function and constraints have to be implemented by deriving from `__interface IFitness`. Besides that, the only class user is supposed to use is `class Pso` which includes all necessary functions for carrying out optimization.
