program exact_solution
implicit none

real*8, parameter :: pii=3.1415926535897932d0

!!variables for elliptic functions and integrals
real*8 :: m            !!elliptic parameter
real*8 :: phi          !!Jacobi amplitude
complex*16 :: u        !!complex argument
complex*16 :: Fc,Ec    !!complete elliptic integrals of first and second kinds
complex*16 :: Fi,Ei    !!incomplete elliptic integrals of first and second kinds 
complex*16 :: sn,cn,dn !!Jacobi elliptic function values

!!variables for physical quantities
real*8 :: x,z,t   !!position and time
real*8 :: lambda  !!aspect ratio
real*8 :: D       !!D value corresponding to the characteristic oribtial
real*8 :: k       !!interface thickness parameter
real*8 :: zI,zR   !!initial z position of interface and reference z for entrainment
real*8 :: RaT,RaC !!thermal and compositional Rayleigh numbers
real*8 :: H       !!internal heating rate

!!numerical resolution variables
integer*4 :: nx,nz    !!mesh size
integer*4 :: nt       !!time series data points
real*8 :: t1,t2       !!initial and final times for time series data

!!arrays
real*8, allocatable :: C_array(:,:),T_array(:,:),H_array(:,:) !!arrays for C, T, and H fields

!!external functions
real*8 :: C,T_func !!functions for C and T
real*8 :: z0       !!function for initial z position of a fluid parcel
real*8 :: vRMS     !!function for RMS velocity

!!testing
character*256 :: fname !!output file name
complex*16 :: kcomplex,InverseJacobiAM
real*8 :: tstart,tfinish
real*8 :: H_horizontal_boundaries
real*8 :: arccot

!!Input Parameters -- note that functions in input_functions.f90 must also be specified
!lambda=1.0d0; k=35.d0; zI=0.5d0; RaT=1.d5; RaC=0.5d5 !!case 1 -- physical parameters
!nx=401; nz=401   !!case 1 -- mesh size
!t1=0.d0; t2=0.01d0 !!case 1 -- time range for entrainment time series
!nt=11             !!case 1 -- # of data points in the entrainment time series

lambda=1.5d0; k=35.d0; zI=0.2d0; RaT=1.d6; RaC=8.d5 !!case 2 -- physical parameters
nx=751; nz=501    !!case 2 -- mesh size
t1=0.d0; t2=0.1d0 !!case 2 -- time range for entrainment time series
nt=11             !!case 2 -- # of data points in the entrainment time series
!!End Input Parameters

!!!complete elliptic integrals
m=1.d9
call complete_elliptic_integrals(m,Fc,Ec)
write(*,*) "Complete Elliptic Integrals:"
write(*,*) "m=",m
write(*,*) "Fc,Ec=",Fc,Ec
write(*,*) ""
!!!end complete elliptic integrals

!!!!!!!!!!!!!!!Incomplete elliptic integrals
phi=1.57d0; m=10.d0
call incomplete_elliptic_integrals(phi,m,Fi,Ei)
write(*,*) "Incomplete Elliptic Integrals:"
write(*,*) "phi,m=",phi,m
write(*,*) "F,E=",Fi,Ei
write(*,*) ""
!!!!!!!!!!!!!!!End Incomplete elliptic integrals

!!!!!!!!!!!!!!!!!!!!!!!!!!!Jacobi elliptic functions
u=1.d0; m=1.d9
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
write(*,*) "Jacobi elliptic functions:"
write(*,*) "u,m=",u,m
write(*,*) "sn,cn,dn=",sn,cn,dn
write(*,*) ""
!!!!!!!!!!!!!!!!!!!!!!!!!!!End Jacobi elliptic functions

!!Functions for Mantle Convection
x=0.2d0; z=0.2d0; t=0.0d0 !!sample coordinate and time for function evaluations
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
!!End Functions for Mantle Convection

!!!!!!Compute H, T, and C arrays and print to file
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
!!!!!!End Compute H, T, and C arrays and print to file

!!!!Entrainment calculation
zR=zI !!reference height
fname="entrainment.dat"
call cpu_time(tstart)
call compute_entrainment(t1,t2,nt,zR,lambda,k,zI,RaT,RaC,nx,nz,fname)
call cpu_time(tfinish)
write(*,*) "Entrainment compute time=",tfinish-tstart
write(*,*) " "
!!!!end Entrainment calculations

!!developer testing
!x=0.1d0; z=0.7d0; t=0.002d0
!call compute_H_func(x,z,t,lambda,k,zI,RaT,RaC,H)
!write(*,*) "Testing compute_H_func:"
!write(*,*) "x,z,t=",x,z,t
!write(*,*) "H=",H

end program exact_solution

