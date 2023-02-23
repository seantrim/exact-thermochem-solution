#example.py

##How to Use
#1) f(t), df/dt, and the integral of f must be specified in input_functions.f90
#2) Run create_library script
#3) Run example.py

import flib ##import library created from f2py3

x=0.999; z=0.999; t=0.005; RaT=1.0e5; RaC=0.5e5; lam=1.0; k=35.0; zI=0.5 ##sample parameters for "temporally periodic" sample problem

##all functions in the Fortran source are available via the flib library -- here are examples for C, T, and H.
##note that Fortran functions are lower case (which may differ from their appearance in the Fortran source files)
##see exact_solution_main.f90 for more functions
Cval=flib.c(x,z,t,lam,k,zI)
Tval=flib.t_func(x,z,t,lam,k,zI,RaT,RaC)
Hval=flib.h_python(x,z,t,lam,k,zI,RaT,RaC)

print("C=",Cval)
print("T=",Tval)
print("H=",Hval)

