# OG
# gfortran -flto -O3 exact_solution_main.f90 exact_solution_routines.f90 ellipFor_modules/elliptic.f90 H_helper_routines.f90 ellipFor_modules/kind_parameters.f90 ellipFor_modules/xelbdj2_all_routines.f90 ellipFor_modules/xgscd_routines.f90 H_func.f90 input_functions.f90 -o exact_solution_code ##maximize speed

# define directory containing ellipFor modules
mod_dir := ellipFor_modules

# GNU Fortran
gfortran: FC = gfortran
gfortran: FOPTS = -flto -O3
gfortran: exact_solution_main
gfortran:
	echo $(FC) $(FOPTS)

# GNU Fortran
gfortran_debug: FC = gfortran
gfortran_debug: FOPTS = -Og -ffpe-trap=invalid,zero,overflow,underflow -ggdb3 -fbounds-check -fbacktrace
gfortran_debug: exact_solution_main
gfortran_debug:
	echo $(FC) $(FOPTS)

exact_solution_main: exact_solution_main.o exact_solution_routines.o exact_solution_routines.mod H_function.o H_function.mod H_helper_routines.o H_helper_routines.mod input_functions.o input_functions.mod elliptic.o elliptic.mod xgscd_routines.o xgscd_routines.mod xelbdj2_all_routines.o xelbdj2_all_routines.mod kind_parameters.o kind_parameters.mod
	$(FC) -o exact_solution_code $(FOPTS) exact_solution_main.f90 exact_solution_routines.f90 H_function.f90 H_helper_routines.f90 input_functions.f90 $(mod_dir)/elliptic.f90 $(mod_dir)/xgscd_routines.f90 $(mod_dir)/xelbdj2_all_routines.f90 $(mod_dir)/kind_parameters.f90

exact_solution_main.o: exact_solution_main.f90 H_function.mod exact_solution_routines.mod kind_parameters.mod
	$(FC) -c $(FOPTS) exact_solution_main.f90 

exact_solution_routines.o exact_solution_routines.mod: exact_solution_routines.f90 input_functions.mod elliptic.mod H_helper_routines.mod H_function.mod kind_parameters.mod
	$(FC) -c $(FOPTS) exact_solution_routines.f90 

H_function.o H_function.mod: H_function.f90 H_helper_routines.mod
	$(FC) -c $(FOPTS) H_function.f90 

H_helper_routines.o H_helper_routines.mod: H_helper_routines.f90 input_functions.mod elliptic.mod
	$(FC) -c $(FOPTS) H_helper_routines.f90 

input_functions.o input_functions.mod: input_functions.f90 kind_parameters.mod
	$(FC) -c $(FOPTS) input_functions.f90 

elliptic.o elliptic.mod: $(mod_dir)/elliptic.f90 kind_parameters.mod xelbdj2_all_routines.mod xgscd_routines.mod
	$(FC) -c $(FOPTS) $(mod_dir)/elliptic.f90

xgscd_routines.o xgscd_routines.mod: $(mod_dir)/xgscd_routines.f90 kind_parameters.mod
	$(FC) -c $(FOPTS) $(mod_dir)/xgscd_routines.f90

xelbdj2_all_routines.o xelbdj2_all_routines.mod: $(mod_dir)/xelbdj2_all_routines.f90 kind_parameters.mod
	$(FC) -c $(FOPTS) $(mod_dir)/xelbdj2_all_routines.f90

kind_parameters.o kind_parameters.mod: $(mod_dir)/kind_parameters.f90
	$(FC) -c $(FOPTS) $(mod_dir)/kind_parameters.f90

clean:
	rm -f *.mod *.o

