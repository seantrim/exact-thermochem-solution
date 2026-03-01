#export FC="ifx" # source /opt/intel/oneapi/setvars.sh
#export FC="gfortran"
#export FFLAGS='-flto -O3' # production

rm -rf builddir
meson setup builddir -Dbuildtype=plain -Dwarning_level=0
meson compile -C builddir -v

