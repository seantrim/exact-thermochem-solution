# Input Options
FC="ifx"     # Fortran compiler: "gfortran" or "ifx" (may need to run source /opt/intel/oneapi/setvars.sh for access to ifx via Intel oneAPI)
mode="production" # build mode: "debug" or "production"

# Automated Build Operations Below
export FC

rm -rf builddir # ensures a new build from scratch

# set up builddir for the selected build mode
if [[ "$mode" == "debug" ]]; then # additional debug options that are not easily selected with built-in meson options
  if [[ "$FC" == "gfortran" ]]; then
    export FFLAGS='-ffpe-trap=invalid,zero,overflow,underflow -fbounds-check -fbacktrace' # floating-point traps and checks for arrays bounds
  elif [[ "$FC" == "ifx" ]]; then
    export FFLAGS='-fpe:0 -check all -fp-model strict' # floating-point traps and checks for arrays bounds and shape violations
  else
    export FFLAGS=''
  fi
  meson setup builddir -Dbuildtype=debug -Dwarning_level=everything # debug
else # apply O3 optimzation to account for meson not applying it at the link stage (bug info: https://github.com/mesonbuild/meson/issues/11318)
  if [[ "$FC" == "gfortran" ]]; then
    export FFLAGS='-O3'
  elif [[ "$FC" == "ifx" ]]; then
    #export FFLAGS='-O3 -ipo' # apply -ipo option because meson does not apply it automatically with built-in link time optimization option 
    export FFLAGS='-O3 -fp-model strict' # neglecting -ipo for now because it is currently not functioning (as of meson 1.7.0 and ifx 2025.3.2 20260112)
  else
    export FFLAGS=''
  fi
  meson setup builddir -Dbuildtype=release -Dwarning_level=0 -Db_lto=true # production
fi

meson compile -C builddir -v       # build executable (artifacts stored in builddir)
mv builddir/exact_solution_code ./ # move executable to current directory
