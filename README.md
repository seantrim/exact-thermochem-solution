# exact-thermochem-solution

## Background

This repository contains files and data supporting the article "Manufacturing an exact solution for 2D thermochemical mantle convection models."

## File Description

Variable Definitions
* $t$ = time
* $C$ = composition
* $T$ = temperature
* $H$ = internal heating rate
* $v_{RMS}$ = root-mean-square velocity
* $E$ = entrainment

SageMath script that can be used to:
* reproduce the formulas for $T$ and $v_{RMS}$
* generate plot data for $C$ and $T$

Maple script for symbolic calculation of $H$

MATLAB script for generating data used in $H$ and $E$

Data files for $E$ corresponding to the "Sample Results" section

Files containing formulas for $H$ (6 files total)
* Fortran
* C++
* Python
* general $f(t)$ (time dependence in the stream function)
* special case of $f(t)=a \sin(\pi b t)$

## Legal

Maple is a trademark of Waterloo Maple Inc.
