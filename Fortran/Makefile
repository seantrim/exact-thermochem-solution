# default compiler
FC := gfortran

# define directory containing ellipFor modules
mod_dir := ellipFor_modules

# GNU Fortran
gfortran: FC := gfortran
gfortran: FOPTS := -flto -O3
gfortran: exact_solution_main

# GNU Fortran (Debug)
gfortran_debug: FC := gfortran
gfortran_debug: FOPTS := -Og -ffpe-trap=invalid,zero,overflow,underflow -Wall -fbounds-check -fbacktrace
gfortran_debug: exact_solution_main

# Intel Fortran
# Note: Intel oneAPI users may have to run the setvars script for access to ifx (e.g., source /opt/intel/oneapi/setvars.sh)
ifx: FC := ifx
ifx: FOPTS := -ipo -O3 -fp-model strict
ifx: exact_solution_main

# Intel Fortran (Debug)
# Note: Intel oneAPI users may have to run the setvars script for access to ifx (e.g., source /opt/intel/oneapi/setvars.sh)
ifx_debug: FC := ifx
ifx_debug: FOPTS := -O0 -g -traceback -fpe:0 -check all -fp-model strict
ifx_debug: exact_solution_main

# Note 1: module names must be lowercase (Fortran compilers produce lowercase module names by default)
# Note 2: object names match the case of the file name

exact_solution_main: exact_solution_main.o exact_solution_routines.o H_function.o H_helper_routines.o input_functions.o elliptic.o xgscd_routines.o xelbdj2_all_routines.o kind_parameters.o
	$(FC) -o exact_solution_code $(FOPTS) exact_solution_main.o exact_solution_routines.o H_function.o H_helper_routines.o input_functions.o elliptic.o xgscd_routines.o xelbdj2_all_routines.o kind_parameters.o

exact_solution_main.o: h_function.mod exact_solution_routines.mod kind_parameters.mod
	$(FC) -c $(FOPTS) exact_solution_main.f90 

exact_solution_routines.o exact_solution_routines.mod: input_functions.mod elliptic.mod h_helper_routines.mod h_function.mod kind_parameters.mod
	$(FC) -c $(FOPTS) exact_solution_routines.f90 

H_function.o h_function.mod: h_helper_routines.mod kind_parameters.mod
	$(FC) -c $(FOPTS) H_function.f90 

H_helper_routines.o h_helper_routines.mod: input_functions.mod elliptic.mod kind_parameters.mod
	$(FC) -c $(FOPTS) H_helper_routines.f90 

input_functions.o input_functions.mod: kind_parameters.mod
	$(FC) -c $(FOPTS) input_functions.f90 

elliptic.o elliptic.mod: kind_parameters.mod xelbdj2_all_routines.mod xgscd_routines.mod
	$(FC) -c $(FOPTS) $(mod_dir)/elliptic.f90

xgscd_routines.o xgscd_routines.mod: kind_parameters.mod
	$(FC) -c $(FOPTS) $(mod_dir)/xgscd_routines.f90

xelbdj2_all_routines.o xelbdj2_all_routines.mod: kind_parameters.mod
	$(FC) -c $(FOPTS) $(mod_dir)/xelbdj2_all_routines.f90

kind_parameters.o kind_parameters.mod:
	$(FC) -c $(FOPTS) $(mod_dir)/kind_parameters.f90

clean:
	rm -f *.mod *.o
