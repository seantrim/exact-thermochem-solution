# exact-thermochem-solution

## Contents
1. [Background](#background)
2. [Variable Definitions](#variable-definitions)
3. [File Description](#file-description)
    1. [SageMath](#sagemath)
    2. [Maple](#maple)
    3. [Fortran](#fortran)
    4. [Data](#data)
4. [How to Use the Fortran Routines](#how-to-use-the-fortran-routines)
    1. [Standalone](#standalone)
    2. [With a Convection Code](#with-a-convection-code)
5. [Legal](#legal)

## Background

This repository contains files and data supporting the article "Manufacturing an exact solution for 2D thermochemical mantle convection models" by S.J. Trim, S.L. Butler, S.S.C. McAdam, and R.J. Spiteri. Computer algebra scripts for the exact solution are provided in SageMath and Maple. Symbolic computation of the internal heating rate is performed using Maple, which has been translated into Fortran. The Fortran routines can be used to calculate quantities from the exact solution both independently and within an existing convection code.

## Variable Definitions
* $t$ = time
* $C$ = composition
* $T$ = temperature
* $H$ = internal heating rate
* $v_{RMS}$ = root-mean-square velocity
* $E$ = entrainment
* $f(t)$ = time dependence of the stream function

## File Description

### [SageMath](/SageMath)
[exact_solution.sage](/SageMath/exact_solution.sage)
* SageMath script that can be used to:
    * symbolically compute the formulas for $T$ and $v_{RMS}$
    * generate data and plots for $C$ and $T$
    * note: does not calculate $H$ (see [Maple](/Maple))

### [Maple](/Maple)
[maple_analytic_solution_include_exterior.mw](/Maple/maple_analytic_solution_include_exterior.mw)
* Maple worksheet for symbolic computation of $H$ (everywhere except at the domain boundaries)
* Results stored in [foo_exterior.m](/Maple/foo_exterior.m) 

[foo_exterior.m](/Maple/foo_exterior.m)
* Results from running [maple_analytic_solution_include_exterior.mw](/Maple/maple_analytic_solution_include_exterior.mw)
* Once loaded into a Maple worksheet, functions such as $C(x,z,t)$, $T(x,z,t)$, and $H(x,a,t)$ become available
    * Saves time compared to rerunning [maple_analytic_solution_include_exterior.mw](/Maple/maple_analytic_solution_include_exterior.mw) 

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

## How to Use the Fortran Routines

### Standalone
Can be used to generate data for $C$, $T$, $H$, $v_{RMS}$, and $E$ from the exact solution.

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
Can be used to calculate $H(x,z,t)$ from within a convection code. These instructions presume that the convection code is written in Fortran. However, options for other programming languages are under development. The following steps are guidelines only. The precise procedure may depend on the particular convection code used.

1. Specify $f(t)$, $df/dt$, $d^2 f/dt^2$ and $\int f dt$ in [input_functions.f90](/Fortran/input_functions.f90)
    * Both cases from the "Sample Results" section are shown as examples in [input_functions.f90](/Fortran/input_functions.f90)
2. Insert calls to the subroutine `compute_H_func`, which returns the value of $H(x,z,t)$, within the source of the convection code where necessary
    * It is presumed that the convection code can accept an internal heating rate the varies in space and time
    * An example of how to call `compute_H_func` is shown in [exact_solution_main.f90](/Fortran/exact_solution_main.f90)
        * The H value is returned in the rightmost argument  
4. Link all f90 files from the [Fortran](/Fortran) folder except [exact_solution_main.f90](/Fortran/exact_solution_main.f90) to the source for the convection code
    * Example: `gfortran -flto -O3 convection_code_source.f90 exact_solution_routines.f90 elliptic.f90 H_helper_routines.f90 xelbdj2_all_routines.f90 xgscd_routines.f90 H_func.f90 input_functions.f90 -o convection_code`
    * In the above example, the source for the convection code is `convection_code_source.f90` and the resulting executable is `convection_code`
        * Modify these names as needed
    * gfortran 11.3.0 or later is recommended
    * Other compilers (and compiler options) may be possible but results should be tested
    * Warning: Duplicate variable/routine names may occur
        * Resolve any related compiler errors
        * Verify that the arguments of `compute_H_func` correspond to the correct values and data types 
5. Run the convection code execuatable as usual

## Legal

This repository is subject to the GPLv3 [license](/LICENSE).

The routines contained within xelbdj2_all_routines.f90 and xgscd_routines.f90 are adapted from routines by Toshio Fukushima available under the CC BY-SA 4.0 license. Original versions of these routines can be found at http://dx.doi.org/10.13140/RG.2.2.27011.66085 and https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum.

Maple is a trademark of Waterloo Maple Inc.
