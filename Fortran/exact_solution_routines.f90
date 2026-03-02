module exact_solution_routines
 use kind_parameters,only: isp,dp
 implicit none
 private
 public :: C,T_func,vRMS,compute_z0,compute_array,create_datafile,compute_entrainment

contains

 real(dp) function vRMS(lambda,t)
  use input_functions,only: f_func
  implicit none
  
  !!input
  real(dp),intent(in) :: lambda !!aspect ratio
  real(dp),intent(in) :: t !!time
  
  !!internal variables
  real(dp), parameter :: pii=3.1415926535897932_dp
  vRMS=pii*sqrt(lambda**2+1._dp)/2._dp/lambda*abs(f_func(t))
 end function vRMS
 
 subroutine compute_entrainment(tmin,tmax,nt,zR,lambda,k,zI,RaT,RaC,nx,nz,fname)
  !!compute entrainment: output is in file named using the input variable fname
  implicit none
  
  !!input
  real(dp)      ,intent(in) :: tmin,tmax !!lower and upper limit for time
  integer(isp)  ,intent(in) :: nt !!number of data points in the time series
  real(dp)      ,intent(in) :: zR !!reference z value
  real(dp)      ,intent(in) :: lambda !!aspect ratio
  real(dp)      ,intent(in) :: k,zI
  real(dp)      ,intent(in) :: RaT,RaC !!Rayleigh numbers
  integer(isp)  ,intent(in) :: nx,nz !!# of mesh points in the x and z directions
  character(256),intent(in) :: fname !!file name
  
  !!internal variables
  real(dp) :: C_array(1:nx,1:nz)
  real(dp) :: t,dt
  real(dp) :: integral,entrainment
  integer(isp) :: kt
  
  open(unit=666,file=fname)
  
  dt=(tmax-tmin)/real(nt-1,dp)
  
  do kt=1,nt
   t=max(dt*real(kt-1,dp),0._dp) !!negative t values are not allowed
   call compute_array('C',t,lambda,k,zI,RaT,RaC,nx,nz,C_array)
   call volume_integral(0._dp,lambda,zR,1._dp,lambda,nx,nz,C_array,integral)
   entrainment=integral/lambda/zI
   write(666,'(2(g16.9))') t,entrainment
  end do
  close(unit=666)
 end subroutine compute_entrainment
 
 subroutine volume_integral(xmin,xmax,zmin,zmax,lambda,nx,nz,array,integral)
  !!computes volume integral of values stored in array over x-z space
  !!Limits of integration given by xmin, xmax, zmin, and zmax.
  implicit none
  
  !!input
  real(dp)    ,intent(in) :: xmin,xmax,zmin,zmax !!limits of integration
  real(dp)    ,intent(in) :: lambda !!aspect ratio
  integer(isp),intent(in) :: nx,nz !!# of mesh points in the x and z directions
  real(dp)    ,intent(in) :: array(1:nx,1:nz) !!array containing the data
  
  !!output
  real(dp),intent(out) :: integral
  
  !!internal variables
  integer(isp) :: kint
  integer(isp) :: kint_min,kint_max
  real(dp)     :: dz
  real(dp)     :: integral_temp(1:nz)
  
  dz=1._dp/real(nz-1,dp)
  
  !integral=0._dp
  kint_min=nint(zmin/dz,isp)+1
  kint_max=nint(zmax/dz,isp)+1
  
  do kint=kint_min,kint_max
   call integral_1D(xmin,xmax,lambda,nx,array(1:nx,kint),integral_temp(kint)) !!create array of integrals over x
  end do
  
  call integral_1D(zmin,zmax,1._dp,nz,integral_temp,integral)
 end subroutine volume_integral
 
 subroutine integral_1D(xmin,xmax,length,nx,array_1D,integral)
  !!compute integral of values stored in array_1D over a single spatial dimension over the interval [0,xmax]
  implicit none
  
  !!input
  real(dp)    ,intent(in) :: xmin,xmax !!lower and upper limits of integration
  real(dp)    ,intent(in) :: length !!length of spatial axis
  integer(isp),intent(in) :: nx !!# of mesh points
  real(dp)    ,intent(in) :: array_1D(1:nx)
  
  !!output
  real(dp),intent(out) :: integral
  
  !!internal variables
  integer(isp) :: iint,iint_min,iint_max
  real(dp)     :: dx
  
  dx=length/real(nx-1,dp)
  
  integral=0._dp
  iint_min=nint(xmin/dx,isp)+1
  iint_max=nint(xmax/dx,isp)+1
  do iint=iint_min,iint_max-1
   integral=integral+dx*(array_1D(iint)+array_1D(iint+1))/2._dp
  end do
 end subroutine integral_1D
 
 subroutine create_datafile(lambda,nx,nz,array,fname)
  !!create datafile for array in column format: x z array
  !!note: output is stored in file named fname
  implicit none
  
  !!input
  real(dp)       ,intent(in) :: lambda !!aspect ratio
  integer(isp)    ,intent(in) :: nx,nz !!# of mesh points in the x and z directions
  real(dp)       ,intent(in) :: array(1:nx,1:nz) !!array containing the data
  character(256),intent(in) :: fname !!output file name
  
  !!internal variables
  integer(isp) :: iint,kint
  real(dp)    :: dx,dz
  real(dp)    :: x,z
  
  open(unit=666,file=fname)
  
  dx=lambda/real(nx-1,dp); dz=1._dp/real(nz-1,dp)
  
  do kint=1,nz
   z=max(min(dz*real(kint-1,dp),1._dp),0._dp)
   do iint=1,nx
    x=max(min(dx*real(iint-1,dp),lambda),0._dp)
    write(666,'(3(g17.9))') x,z,array(iint,kint)
   end do
  end do
  
  close(666)
 end subroutine create_datafile
 
 subroutine compute_array(option,t,lambda,k,zI,RaT,RaC,nx,nz,array)
  !!compute array of function values for a time t
  use H_function,only: compute_H_func
  implicit none
  
  !!input
  character(1),intent(in) :: option !!specify function to be computed -- valid options are: "C", "T", "D", and "H".
  real(dp)     ,intent(in) :: t !!time
  real(dp)     ,intent(in) :: lambda,k,zI !!aspect ratio, sharpness parameter, and z value of interface at t=0
  real(dp)     ,intent(in) :: RaT,RaC !!Rayleigh numbers
  integer(isp)  ,intent(in) :: nx,nz !!# of mesh points in the x and z directions
  
  !!output
  real(dp),intent(out) :: array(1:nx,1:nz)
  
  !!internal variables
  integer(isp) :: iint,kint
  real(dp) :: dx,dz
  real(dp) :: x,z
  real(dp) :: H_func
  
  dx=lambda/real(nx-1,dp); dz=1._dp/real(nz-1,dp)
  
  if (option.eq.'C') then
   do kint=1,nz
    z=max(min(dz*real(kint-1,dp),1._dp),0._dp)
    do iint=1,nx
     x=max(min(dx*real(iint-1,dp),lambda),0._dp)
     array(iint,kint)=C(x,z,t,lambda,k,zI)
    end do
   end do
  else if (option.eq.'T') then
   do kint=1,nz
    z=max(min(dz*real(kint-1,dp),1._dp),0._dp)
    do iint=1,nx
     x=max(min(dx*real(iint-1,dp),lambda),0._dp)
     array(iint,kint)=T_func(x,z,t,lambda,k,zI,RaT,RaC)
    end do
   end do
  else if (option.eq.'D') then
   do kint=1,nz
    z=max(min(dz*real(kint-1,dp),1._dp),0._dp)
    do iint=1,nx
     x=max(min(dx*real(iint-1,dp),lambda),0._dp)
     array(iint,kint)=D(x,z,lambda)
    end do
   end do
  else if (option.eq.'H') then
   do kint=1,nz
    z=max(min(dz*real(kint-1,dp),1._dp),0._dp)
    do iint=1,nx
     x=max(min(dx*real(iint-1,dp),lambda),0._dp)
     call compute_H_func(x,z,t,lambda,k,zI,RaT,RaC,H_func)
     array(iint,kint)=H_func
    end do
   end do
  else
   write(*,*) "Error in compute_array: unrecognized option."
  end if
 end subroutine compute_array
 
 real(dp) function T_func(x,z,t,lambda,k,zI,RaT,RaC)
  !!temperature function
  !!assumes x belongs to (-lambda/2,3/2*lambda) and z belongs to (-1,2)
  use input_functions,only: f_func
  implicit none
  
  !!inputs
  real(dp),intent(in) :: x,z,t  !!position and time
  real(dp),intent(in) :: lambda !!aspect ratio
  real(dp),intent(in) :: k,zI   !!sharpness parameter and initial layer height
  real(dp),intent(in) :: RaT,RaC !!Rayleigh numbers
  
  !!internal variables
  real(dp), parameter :: pii=3.1415926535897932_dp
  
  T_func=(-(pii**3*(lambda**2+1._dp)**2/lambda**3)*cos(pii*x/lambda)*sin(pii*z)*f_func(t)+RaC*C(x,z,t,lambda,k,zI)+&
  &(RaT-RaC)*(1._dp-z))/RaT
 end function T_func
 
 real(dp) function C(x,z,t,lambda,k,zI)
  !!compute composition
  !!assumes x belongs to (-lambda/2,3/2*lambda) and z belongs to (-1,2)
  implicit none
  
  !!inputs
  real(dp),intent(in) :: x,z,t  !!position and time
  real(dp),intent(in) :: lambda !!aspect ratio
  real(dp),intent(in) :: k,zI   !!sharpness parameter and initial layer height
  
  !!internal variables
  real(dp) :: z0
  
  call compute_z0(x,z,t,lambda,z0)
  C=1._dp/(1._dp+exp(-2._dp*k*(zI-z0)))
 end function C
 
 subroutine compute_z0(x,z,t,lambda,z0)
  !!Compute z0: the initial z value for a fluid parcel at position (x,z) at time t.
  !!assumes x belongs to (-lambda/2,3/2*lambda) and z belongs to (-1,2)
  use input_functions,only: f_integral
  use elliptic,only: incomplete_elliptic_integrals,Jacobi_elliptic_functions
  use H_helper_routines,only: arccot
  implicit none
  
  !!inputs
  real(dp),intent(in) :: x,z,t  !!position and time
  real(dp),intent(in) :: lambda !!aspect ratio
  
  !!output
  real(dp),intent(out) :: z0
  
  !!internal variables
  real(dp), parameter :: pii=3.1415926535897932_dp
  real(dp) :: Q,arg,eQ,bigZ0 !!variables for sidewalls
  real(dp) :: phi,m !!elliptic amplitude and parameter
  real(dp) :: x_,z_
  complex(dp) :: F,E !!elliptic integrals of the first and second kinds
  complex(dp) :: u   !!argument for the Jacobi elliptic functions
  complex(dp) :: sn,cn,dn !!Jacobi elliptic function values
  
  if ((z.eq.0._dp).or.(z.eq.1._dp).or.((x.eq.(lambda/2._dp)).and.(z.eq.0.5_dp))) then
   z0=z
  else if ((x.eq.0._dp).or.(x.eq.lambda)) then
   arg=pii*z
   Q=log(abs(1._dp/sin(arg)+1._dp/tan(arg)))-pii/lambda*S(x,lambda)*f_integral(t)
   eQ=exp(Q)
   bigZ0=(eQ-1._dp/eQ)/2._dp
   if (bigZ0.ge.0._dp) then
    z0=arccot(bigZ0)/pii
   else
    z0=1._dp+arccot(bigZ0)/pii
   end if
  else
   if (x.lt.0._dp) then !!beyond left sidewall
    x_=-x
   else if (x.gt.lambda) then !!beyond right sidewall
    x_=2._dp*lambda-x
   else !!domain interior
    x_=x
   end if
   if (z.lt.0._dp) then !!beyond bottom
    z_=-z
   else if (z.gt.1._dp) then !!beyond top
    z_=2._dp-z  
   else !!domain interior
    z_=z
   end if
   phi=pii*z_; m=1._dp/D(x_,z_,lambda)**2._dp
   call incomplete_elliptic_integrals(phi,m,F,E)
   u=F
   if (t.ne.0._dp) u%IM=u%IM-S(x_,lambda)*(pii**2._dp*D(x_,z_,lambda)/lambda)*f_integral(t)
   call Jacobi_elliptic_functions(u,m,sn,cn,dn)
   if (abs(cn%RE).gt.1._dp) then
    write(*,*) "compute_z0: large cn magnitude detected -- cn=",cn%RE
    stop
   end if
   z0=acos(cn%RE)/pii
  
  
  ! phi=pii*z; m=1._dp/D(x,z,lambda)**2._dp
  ! call incomplete_elliptic_integrals(phi,m,F,E)
  ! u=F
  ! if (t.ne.0._dp) u%IM=u%IM-S(x,lambda)*(pii**2._dp*D(x,z,lambda)/lambda)*f_integral(t)
  ! call Jacobi_elliptic_functions(u,m,sn,cn,dn)
  ! if (abs(cn%RE).gt.1._dp) then
  !  write(*,*) "compute_z0: large cn magnitude detected -- cn=",cn%RE
  !  stop
  ! end if
  ! z0=acos(cn%RE)/pii
  end if
 end subroutine compute_z0
 
 real(dp) function D(x,z,lambda)
  !!characteristic orbital value
  implicit none
  
  !!inputs
  real(dp),intent(in) :: x,z,lambda
  
  !!internal variables
  real(dp), parameter :: pii=3.1415926535897932_dp
  
  D=abs(sin(pii*x/lambda)*sin(pii*z))
 end function D
  
 real(dp) function S(x,lambda)
  !!transformed unit step function
  !!assumes -lambda/2 <= x <= 3/2*lambda
  implicit none
  real(dp),intent(in) :: x,lambda

  if (x.le.lambda/2._dp) then
   S=1._dp
  else
   S=-1._dp
  end if
 end function S

end module exact_solution_routines 
