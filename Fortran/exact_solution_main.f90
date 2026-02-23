program exact_solution
 use kind_parameters,only: isp,dp
 use H_function,only: compute_H_func
 use exact_solution_routines,only: C,T_func,vRMS,compute_z0,compute_array,create_datafile,compute_entrainment
 implicit none
 
 !!Variables for Physical Quantities
 real(dp) :: x,z,t   !!position and time
 real(dp) :: lambda  !!aspect ratio
 real(dp) :: k       !!interface thickness parameter
 real(dp) :: z0      !!initial z position of a fluid parcel
 real(dp) :: zI,zR   !!initial z position of interface and reference z for entrainment
 real(dp) :: RaT,RaC !!thermal and compositional Rayleigh numbers
 real(dp) :: H       !!internal heating rate
 
 !!Numerical Resolution Variables
 integer(isp) :: nx,nz    !!mesh size
 integer(isp) :: nt       !!time series data points
 real(dp)     :: t1,t2    !!initial and final times for time series data
 
 !!Arrays
 real(dp), allocatable :: C_array(:,:),T_array(:,:),H_array(:,:) !!arrays for C, T, and H fields
 
 !!Internal Variables
 character(256) :: fname          !!output file name
 real(dp)       :: tstart,tfinish !!compute time variables
 
 !!!!Input Parameters -- note that functions in input_functions.f90 must also be specified
 lambda=1.0_dp; k=35._dp; zI=0.5_dp; RaT=1.e5_dp; RaC=0.5e5_dp !!case 1 -- physical parameters
 nx=401; nz=401     !!case 1 -- mesh size
 t1=0._dp; t2=0.01_dp !!case 1 -- time range for entrainment time series
 nt=11              !!case 1 -- # of data points in the entrainment time series
 
 !lambda=1.5_dp; k=35._dp; zI=0.2_dp; RaT=1.e6_dp; RaC=8.e5_dp !!case 2 -- physical parameters
 !nx=751; nz=501    !!case 2 -- mesh size
 !t1=0._dp; t2=0.1_dp !!case 2 -- time range for entrainment time series
 !nt=11             !!case 2 -- # of data points in the entrainment time series
 !!End Input Parameters
 
 !!!!Functions for Mantle Convection
 x=0.999_dp; z=0.999_dp; t=0.005_dp !!sample coordinate and time for function evaluations
 call compute_z0(x,z,t,lambda,z0)
 call compute_H_func(x,z,t,lambda,k,zI,RaT,RaC,H)
 print *, "Functions for Mantle Convection:"
 print *, "x,z,t=",x,z,t
 print *, "lambda,k,zI,RaT,RaC=",lambda,k,zI,RaT,RaC
 print *, "vRMS=",vRMS(lambda,t)
 print *, "z0=",z0
 print *, "C=",C(x,z,t,lambda,k,zI)
 print *, "T=",T_func(x,z,t,lambda,k,zI,RaT,RaC)
 print *, "H=",H
 print *, " "
 
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
