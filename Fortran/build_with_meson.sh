#export FC="ifx" # source /opt/intel/oneapi/setvars.sh
export FC="gfortran"
#export FFLAGS='-flto -O3' # production
#export FFLAGS='-Og -ffpe-trap=invalid,zero,overflow,underflow -ggdb3 -fbounds-check -fbacktrace' # debug

rm -rf builddir # ensures a new build from scratch
#meson setup builddir -Dbuildtype=release -Dwarning_level=0 -Db_lto=true # production
meson setup builddir -Dbuildtype=debug -Dwarning_level=everything # debug
meson compile -C builddir -v
mv builddir/exact_solution_code ./
