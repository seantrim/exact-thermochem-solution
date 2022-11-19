# exact-thermochem-solution

## Contents
1. [Background](#background)
2. [Variable Definitions](#variable-definitions)
3. [File Description](#file-description)
    1. [SageMath](#sagemath)
    2. [Maple](#maple)
    3. [Fortran](#fortran)
    4. [Data](#data)
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
* f(t) = time dependence of the stream function

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
* The result was adapted to free form Fortran for use in [H_func.f90](/Fortran/H_func.f90)

### [Fortran](/Fortran)
[exact_solution_main.f90](/Fortran/exact_solution_main.f90)
* Main program for [standalone](#standalone) version of the Fortran routines

[exact_solution_routines.f90](/Fortran/exact_solution_routines.f90)
* Contains routines for calculating physical quantities including:
    * $C$, $T$, $H$, $v_{RMS}$, and $E$

[H_func.f90](/Fortran/H_func.f90)
* Contains routine for evaluating $H(x,z,t)$

[H_helper_routines.f90](/Fortran/H_helper_routines.f90)
* Contains Fortran equivalents of Maple functions referenced in [H_func.f90](/Fortran/H_func.f90)

[input_functions.f90](/Fortran/input_functions.f90)
* Contains user defined function $f(t)$, with its integral and derivatives

[elliptic.f90](/Fortran/elliptic.f90)
* Contains routines that evaluate elliptic integrals and Jacobi elliptic functions for the input parameter ranges needed
* Input ranges were generalized by combining identities with calls to routines from [xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90) and [xgscd_routines.f90](/Fortran/xgscd_routines.f90)

[xelbdj2_all_routines.f90](/Fortran/xelbdj2_all_routines.f90)
* Contains routines for evaluation of associate incomplete elliptic integrals of first, second, and third kinds
* Assumes standard input parameter ranges
* Adapted from routines by Toshio Fukushima (see [Legal](#legal))

[xgscd_routines.f90](/Fortran/xgscd_routines.f90)
* Contains routines for the evaluation of the Jacobi elliptic functions sn, cn, and dn
* Assumes standard input parameter ranges
* Adapted from routines by Toshio Fukushima (see [Legal](#legal))

[compile_script](/Fortran/compile_script)
* Terminal script for compiling the [standalone](#standalone) version of the Fortran routines using gfortran

[exact_solution_code](/Fortran/exact_solution_code)
* Sample exexcutable for the [standalone](#standalone) version of the Fortran routines
* Results from running [compile_script](/Fortran/compile_script) in the terminal using the .f90 files in the [Fortran](/Fortran) folder
* gfortran 11.3.0 was used  

### [Data](/Data)
[entrainment_sample_1_401x401.dat](/Data/entrainment_sample_1_401x401.dat)
* $E$ time series data for temporally periodic case in "Sample Results" section 
* Computed using [Fortran](/Fortran) routines

[entrainment_sample_2_751x501.dat](/Data/entrainment_sample_2_751x501.dat)
* $E$ time series data for approaching steady state case in "Sample Results" section
* Computed using [Fortran](/Fortran) routines

## How to Use Fortran Routines

### Standalone
Can be used to generate data for $C$, $T$, $H$, $v_{RMS}$, and $E$ from the exact solution

1. Specify $f(t)$, $df/dt$, $d^2 f/dt^2$ and $\int f dt$ in [input_functions.f90](/Fortran/input_functions.f90)
    * Both cases from the "Sample Results" section are shown as examples in [input_functions.f90](/Fortran/input_functions.f90)
2. Specify physical parameters (e.g., aspect ratio, Rayleigh numbers, etc.) in [exact_solution_main.f90](/Fortran/exact_solution_main.f90)
    * Search for "!!Input Parameters" in the comments
    * Both cases from the "Sample Results" section are shown as examples
3. Compile the code by running [compile_script](/Fortran/compile_script) from the command line.
    * Linux: `$ source compile_script`
        * This produces the executable [exact_solution_code](/Fortran/exact_solution_code)
    * gfortran 11.3.0 or later is recommended
    * Other compilers may be possible but results should be tested
4. Run [exact_solution_code](/Fortran/exact_solution_code) from the command line
    * Linux: `$ ./exact_solution_code`
    * This will produce data for snapshots of $C$, $T$, and $H$
    * This will also produce data for a time series of $E$
    * Example routine calls for several quantities are shown in [exact_solution_main.f90](/Fortran/exact_solution_main.f90)
    * Modifying the examples in [exact_solution_main.f90](/Fortran/exact_solution_main.f90) can be done to customize the output

### With a Convection Code
* Calculate $H(x,z,t)$ to be used in the advection--diffusion equation for $T$. 

## Legal

This repository is subject to the GPLv3 [license](/LICENSE).

The routines contained within xelbdj2_all_routines.f90 and xgscd_routines.f90 are adapted from routines by Toshio Fukushima available under the CC BY-SA 4.0 license. Original versions of these routines can be found at http://dx.doi.org/10.13140/RG.2.2.27011.66085 and https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum.

Maple is a trademark of Waterloo Maple Inc.
