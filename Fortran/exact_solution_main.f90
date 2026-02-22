program exact_solution
 use H_function,only: compute_H_func
 use exact_solution_routines,only: C,T_func,vRMS,compute_z0,compute_array,create_datafile,compute_entrainment
 implicit none
 
 !!Variables for Physical Quantities
 real*8 :: x,z,t   !!position and time
 real*8 :: lambda  !!aspect ratio
 real*8 :: k       !!interface thickness parameter
 real*8 :: z0      !!initial z position of a fluid parcel
 real*8 :: zI,zR   !!initial z position of interface and reference z for entrainment
 real*8 :: RaT,RaC !!thermal and compositional Rayleigh numbers
 real*8 :: H       !!internal heating rate
 
 !!Numerical Resolution Variables
 integer*4 :: nx,nz    !!mesh size
 integer*4 :: nt       !!time series data points
 real*8 :: t1,t2       !!initial and final times for time series data
 
 !!Arrays
 real*8, allocatable :: C_array(:,:),T_array(:,:),H_array(:,:) !!arrays for C, T, and H fields
 
 !!Internal Variables
 character*256 :: fname   !!output file name
 real*8 :: tstart,tfinish !!compute time variables
 
 !!!!Input Parameters -- note that functions in input_functions.f90 must also be specified
 lambda=1.0d0; k=35.d0; zI=0.5d0; RaT=1.d5; RaC=0.5d5 !!case 1 -- physical parameters
 nx=401; nz=401     !!case 1 -- mesh size
 t1=0.d0; t2=0.01d0 !!case 1 -- time range for entrainment time series
 nt=11              !!case 1 -- # of data points in the entrainment time series
 
 !lambda=1.5d0; k=35.d0; zI=0.2d0; RaT=1.d6; RaC=8.d5 !!case 2 -- physical parameters
 !nx=751; nz=501    !!case 2 -- mesh size
 !t1=0.d0; t2=0.1d0 !!case 2 -- time range for entrainment time series
 !nt=11             !!case 2 -- # of data points in the entrainment time series
 !!End Input Parameters
 
 !!!!Functions for Mantle Convection
 x=0.999d0; z=0.999d0; t=0.005d0 !!sample coordinate and time for function evaluations
 call compute_z0(x,z,t,lambda,z0)
 call compute_H_func(x,z,t,lambda,k,zI,RaT,RaC,H)
 write(*,*) "Functions for Mantle Convection:"
 write(*,*) "x,z,t=",x,z,t
 write(*,*) "lambda,k,zI,RaT,RaC=",lambda,k,zI,RaT,RaC
 write(*,*) "vRMS=",vRMS(lambda,t)
 write(*,*) "z0=",z0
 write(*,*) "C=",C(x,z,t,lambda,k,zI)
 write(*,*) "T=",T_func(x,z,t,lambda,k,zI,RaT,RaC)
 write(*,*) "H=",H
 write(*,*) " "
 
 !!!!Compute H, T, and C arrays and print to file
 allocate(H_array(1:nx,1:nz),T_array(1:nx,1:nz),C_array(1:nx,1:nz))
 call cpu_time(tstart)
 call compute_array('H',t,lambda,k,zI,RaT,RaC,nx,nz,H_array)
 call cpu_time(tfinish)
 fname="H_data.dat"
 call create_datafile(lambda,nx,nz,H_array,fname)
 call compute_array('T',t,lambda,k,zI,RaT,RaC,nx,nz,T_array)
 fname="T_data.dat"
 call create_datafile(lambda,nx,nz,T_array,fname)
 call compute_array('C',t,lambda,k,zI,RaT,RaC,nx,nz,C_array)
 fname="C_data.dat"
 call create_datafile(lambda,nx,nz,C_array,fname)
 write(*,*) "H array compute time=",tfinish-tstart
 write(*,*) " "
 
 !!!!Entrainment calculation
 zR=zI !!reference height
 fname="entrainment.dat"
 call cpu_time(tstart)
 call compute_entrainment(t1,t2,nt,zR,lambda,k,zI,RaT,RaC,nx,nz,fname)
 call cpu_time(tfinish)
 write(*,*) "Entrainment compute time=",tfinish-tstart
 write(*,*) " "
 
end program exact_solution
