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
complex*16 :: kcomplex,InverseJacobiAM

!!Input Parameters -- note that functions in input_functions.f90 must also be specified
!lambda=1.0d0; k=35.d0; zI=0.5d0; RaT=1.d5; RaC=1.d5 !!case 1 -- physical parameters
!nx=151; nz=151   !!case 1 -- mesh size
t1=0.d0; t2=0.01d0 !!case 1 -- time range for entrainment time series
nt=11             !!case 1 -- # of data points in the entrainment time series

lambda=1.5d0; k=35.d0; zI=0.2d0; RaT=1.d6; RaC=8.d5 !!case 2 -- physical parameters
!nx=601; nz=401    !!case 2 -- mesh size 
nx=301; nz=201    !!case 2 -- mesh size
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
x=0.05d0; z=0.1d00; t=0.00000d0; !!sample coordinate and time for function evaluations
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

!!!!!!H array test -- use parameters from H test above
allocate(H_array(1:nx,1:nz),T_array(1:nx,1:nz),C_array(1:nx,1:nz))
call compute_array('H',t,lambda,k,zI,RaT,RaC,nx,nz,H_array)
call create_datafile(lambda,nx,nz,H_array,"H_data.dat")
call compute_array('T',t,lambda,k,zI,RaT,RaC,nx,nz,T_array)
call create_datafile(lambda,nx,nz,T_array,"T_data.dat")
call compute_array('C',t,lambda,k,zI,RaT,RaC,nx,nz,C_array)
call create_datafile(lambda,nx,nz,C_array,"C_data.dat")
!!!!!!End H array test

!!!!Entrainment calculation
zR=zI !!reference height
call compute_entrainment(t1,t2,nt,zR,lambda,k,zI,RaT,RaC,nx,nz,"entrainment.dat")
!!!!end Entrainment calculations

!!For testing -- 
!write(*,*) "Testing: verify accuracy with Maple"
!u=(3.0892327760299629d0,0.d0); kcomplex=(19.107322609297285d0,0.d0)
!write(*,*) "u,kcomplex=",u,kcomplex
!write(*,*) "InverseJacobiAM(u,kcomplex)=",InverseJacobiAM(u,kcomplex)
!write(*,*) "Maple value = 0.082265508499299800 - 0.45413520599782981 I"
!m=kcomplex%RE**2
!call incomplete_elliptic_integrals(u%RE,m,Fi,Ei)
!write(*,*) "Fi=",Fi
!write(*,*) "Ei=",Ei
!write(*,*) ""

end program exact_solution

