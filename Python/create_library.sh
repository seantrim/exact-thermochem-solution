#script for creating Python library from Fortran routines for H and related functions
f2py3 -c --fcompiler=gfortran --f90flags=-flto --opt=-O3 -m flib ../Fortran/*.f90

