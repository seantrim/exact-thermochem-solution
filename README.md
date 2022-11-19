# exact-thermochem-solution

## Contents
1. [Background](#background)
2. [Variable Definitions](#variable-definitions)
3. [File Description](#file-description)
4. [How to Use Fortran Routines](#how-to-use-fortran-routines)
    1. [Standalone](#standalone)
    2. [With a Convection Code](#with-a-convection-code)
5. [Legal](#legal)

## Background

This repository contains files and data supporting the article "Manufacturing an exact solution for 2D thermochemical mantle convection models" by S.J. Trim, S.L. Butler, S.S.C. McAdam, and R.J. Spiteri.

## Variable Definitions
* $t$ = time
* $C$ = composition
* $T$ = temperature
* $H$ = internal heating rate
* $v_{RMS}$ = root-mean-square velocity
* $E$ = entrainment

## File Description

### [SageMath](/SageMath)
[exact_solution.sage](/SageMath/exact_solution.sage)
* SageMath script that can be used to:
    * symbolically compute the formulas for $T$ and $v_{RMS}$
    * generate plot data for $C$ and $T$

### [Maple](/Maple)
[maple_analytic_solution.mw](/Maple/maple_analytic_solution.mw)
* Maple worksheet for symbolic computation of $H$
* Results stored in [foo.m](/Maple/foo.m) 

[foo.m](/Maple/foo.m)
* Results from running [maple_analytic_solution.mw](/Maple/maple_analytic_solution.mw)
* Once loaded into a Maple worksheet, functions such as $C(x,z,t)$, $T(x,z,t)$, and $H(x,a,t)$ become available
    * Saves time compared to rerunning [maple_analytic_solution.mw](/Maple/maple_analytic_solution.mw) 

[Fortran_code_generation.mw](/Maple/Fortran_code_generation.mw)
* Maple worksheet for translating Maple's $H$ formula into Fortran 77 code
* The result was adapted to free form Fortran for use in H_func.f90

### [Fortran](/Fortran)
* Fortran code that can calculate $C$, $T$, $H$, $v_{RMS}$, and $E$.

### [Data](/Data)
[entrainment_sample_1_401x401.dat](/Data/entrainment_sample_1_401x401.dat)
* $E$ time series data for temporally periodic case in "Sample Results" section 
* Computed using [Fortran](/Fortran) routines

[entrainment_sample_2_751x501.dat](/Data/entrainment_sample_2_751x501.dat)
* $E$ time series data for approaching steady state case in "Sample Results" section
* Computed using [Fortran](/Fortran) routines

## How to Use Fortran Routines

### Standalone
* Calculate $C$, $T$, $H$, $v_{RMS}$, and $E$.

### With a Convection Code
* Calculate $H(x,z,t)$ to be used in the advection--diffusion equation for $T$. 

## Legal

This repository is subject to the GPLv3 [license](/LICENSE).

The routines contained within xelbdj2_all_routines.f90 and xgscd_routines.f90 are adapted from the routines by Toshio Fukushima available under the CC BY-SA 4.0 license. Original versions of these routines can be found at http://dx.doi.org/10.13140/RG.2.2.27011.66085 and https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum.

Maple is a trademark of Waterloo Maple Inc.
