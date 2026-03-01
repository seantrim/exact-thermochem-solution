project('exact_solution_code', 'fortran', default_options : ['buildtype=plain'])
add_project_arguments(language : ['fortran'])
executable('exact_solution_code', files('exact_solution_main.f90', 'exact_solution_routines.f90', 'H_function.f90', 'H_helper_routines.f90', 'input_functions.f90', 'ellipFor_modules/elliptic.f90', 'ellipFor_modules/xelbdj2_all_routines.f90', 'ellipFor_modules/xgscd_routines.f90', 'ellipFor_modules/kind_parameters.f90'))
