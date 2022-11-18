# exact-thermochem-solution
1. [Background](#background)
2. [Variable Definitions](#variable-definitions)
3. [File Description](#file-description)
4. [Legal](#legal)

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

[exact_solution.sage](/exact_solution.sage)
* SageMath script that can be used to:
    * symbolically compute the formulas for $T$ and $v_{RMS}$
    * generate plot data for $C$ and $T$

[Maple](/Maple)
* Maple script for symbolic calculation of $H$

[Fortran](/Fortran)
* Fortran code that can calculate $C$, $T$, $H$, $v_{RMS}$, and $E$.

Data
* Data files for $E$ corresponding to the "Sample Results" section
* Computed using Fortran routines

## Legal

This repository is subject to the GPLv3 [license](/LICENSE).

The routines contained within xelbdj2_all_routines.f90 and xgscd_routines.f90 are adapted from the routines by Toshio Fukushima available under the CC BY-SA 4.0 license. Original versions of these routines can be found at http://dx.doi.org/10.13140/RG.2.2.27011.66085 and https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum.

Maple is a trademark of Waterloo Maple Inc.
