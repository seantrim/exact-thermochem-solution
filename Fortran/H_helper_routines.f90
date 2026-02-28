module H_helper_routines
 use kind_parameters,only: dp
 implicit none
 private

 public :: csc,cot,arccot

 public :: JacobiZeta

 public :: F,dFdt,d2Fdt2

 public :: signum,Heaviside,Dirac,Dirac_derivative

 public :: InverseJacobiAM
 public :: compute_JacobiSN_CN_DN,JacobiSN,JacobiCN,JacobiDN
 public :: EllipticK,EllipticE,compute_EllipticK_EllipticE

contains
 
 real(dp) function csc(theta)
  !!cosecant function
  implicit none
  
  !!input
  real(dp),intent(in) :: theta
  
  csc=1._dp/sin(theta)
 end function csc
 
 real(dp) function cot(theta)
  !!cotangent function 
  implicit none
  
  !!input
  real(dp),intent(in) :: theta
  
  cot=1._dp/tan(theta)
 end function cot

 real(dp) function arccot(x)
  !!inverse cotangent function -- range is (-pi/2,pi/2]-{0} for computational convenience
  !!assuming x is real
  implicit none
   
  !!input
  real(dp),intent(in) :: x
    
  !!internal variables
  real(dp), parameter :: pii=3.1415926535897932_dp
   
  if (x.lt.0._dp) then
   arccot=-pii/2._dp-atan(x)
  else
   arccot=pii/2._dp-atan(x)
  end if
 end function arccot
 
 complex(dp) function JacobiZeta(u,k)
  !!Jacobi Zeta function
  !!assumes imaginary part of k is zero
  implicit none
  
  !!input
  complex(dp),intent(in) :: k !!elliptic modulus
  complex(dp),intent(in) :: u !!argument (equivalent to F(phi|m))
  
  call compute_JacobiZeta(u,k,JacobiZeta)

  contains

   subroutine compute_JacobiZeta(u,k,JacobiZeta)
    !!Jacobi Zeta function
    !!assumes imaginary part of k is zero
    !!assumes the Jacobi amplitude (phi) is real
    use elliptic,only: complete_elliptic_integrals,incomplete_elliptic_integrals,Jacobi_elliptic_functions
    implicit none
    
    !!input
    complex(dp),intent(in) :: k !!elliptic modulus
    complex(dp),intent(in) :: u !!argument (equivalent to F(phi|m))
    
    !!output
    complex(dp),intent(out) :: JacobiZeta
    
    !!internal variables
    real(dp), parameter :: tol_IM=1.e-3_dp !!tolerance for cn%IM size
    real(dp)    :: m !!elliptic parameter
    real(dp)    :: phi !!Jacobi amplitude
    complex(dp) :: Fc,Ec,F,E !!elliptic integral values
    complex(dp) :: twoFc !!two times Fc
    complex(dp) :: sn,cn,dn !!Jacobi elliptic function values
    real(dp)    :: Nm_IM,Np_IM,Nm_RE,Np_RE !!periodicity multiplier
    real(dp)    :: fraction_Nm_RE,fraction_Np_RE,fraction_Nm_IM,fraction_Np_IM
    real(dp)    :: delta_Nm_RE,delta_Np_RE,delta_Nm_IM,delta_Np_IM
    logical     :: IM_periodic,RE_periodic
    logical     :: int_p_RE,int_m_RE,int_p_IM,int_m_IM
    
    m=k%RE**2
    call complete_elliptic_integrals(m,Fc,Ec)
    twoFc=2._dp*Fc !!store this for frequent use
    call Jacobi_elliptic_functions(u,m,sn,cn,dn)
    
    if (abs(cn%RE).gt.1._dp) then
     write(*,*) "Error in compute_JacobiZeta -- abs(cn%RE) is greater than unity"
     write(*,*) "u,k,m=",u,k,m
     write(*,*) "cn=",cn
     stop
    else if (abs(cn%IM).gt.tol_IM) then
     write(*,*) "Error in compute_JacobiZeta -- cn%IM is not zero"
     write(*,*) "u,k,m=",u,k,m
     write(*,*) "cn=",cn
     stop
    end if
    phi=acos(cn%RE) !!produces phi belonging to [0,pi] which works with incomplete_elliptic_integrals routine
    call incomplete_elliptic_integrals(phi,m,F,E) !!assuming imaginary part of phi is negligible
    
    !!Re[u] not in the presumed range? 
    RE_periodic=.not.((((0._dp.le.u%RE).and.(u%RE.le.(twoFc%RE))).or.(((twoFc%RE).lt.u%RE).and.(u%RE.le.0._dp))))
    if (RE_periodic.eqv..true.) then
     Nm_RE=(u%RE-F%RE)/(twoFc%RE)
     Np_RE=(u%RE+F%RE)/(twoFc%RE)
     fraction_Nm_RE=abs(mod(Nm_RE,1._dp)) !fractional parts of N values -- should be very close to zero or unity for integer values of N
     fraction_Np_RE=abs(mod(Np_RE,1._dp))
     delta_Nm_RE=min(fraction_Nm_RE,1._dp-fraction_Nm_RE)
     delta_Np_RE=min(fraction_Np_RE,1._dp-fraction_Np_RE)
     if (delta_Np_RE.lt.delta_Nm_RE) then
      int_p_RE=.true.; int_m_RE=.false.
     else
      int_p_RE=.false.; int_m_RE=.true.
     end if
     if (int_p_RE.eqv..true.) then
      E%RE=2._dp*Ec%RE*Np_RE-E%RE !!use periodicity of elliptic integrals
     else if (int_m_RE.eqv..true.) then
      E%RE=2._dp*Ec%RE*Nm_RE+E%RE
     end if
    end if
    
    !!Im[u] not in the presumed range?
    IM_periodic=.not.((((0._dp.le.u%IM).and.(u%IM.le.(twoFc%IM))).or.(((twoFc%IM).lt.u%IM).and.(u%IM.le.0._dp))))
    if (IM_periodic.eqv..true.) then
     Nm_IM=(u%IM-F%IM)/(twoFc%IM)
     Np_IM=(u%IM+F%IM)/(twoFc%IM)
     fraction_Nm_IM=abs(mod(Nm_IM,1._dp)) !!fractional parts of N values -- should be very close to zero or unity for integer values of N
     fraction_Np_IM=abs(mod(Np_IM,1._dp))
     delta_Nm_IM=min(fraction_Nm_IM,1._dp-fraction_Nm_IM)
     delta_Np_IM=min(fraction_Np_IM,1._dp-fraction_Np_IM)
     if (delta_Np_IM.lt.delta_Nm_IM) then
      int_p_IM=.true.; int_m_IM=.false.
     else
      int_p_IM=.false.; int_m_IM=.true.
     end if
     if (int_p_IM.eqv..true.) then
      E%IM=2._dp*Ec%IM*Np_IM-E%IM
     else if (int_m_IM.eqv..true.) then
      E%IM=2._dp*Ec%IM*Nm_IM+E%IM
     end if
    end if
    
    JacobiZeta=E-Ec*u/Fc
   end subroutine compute_JacobiZeta

 end function JacobiZeta
 
 real(dp) function F(t)
  !!wrapper function for integral of f(t)
  use input_functions,only: f_integral
  implicit none
  real(dp),intent(in) :: t
  F=f_integral(t)
 end function F
 
 real(dp) function dFdt(t)
  !!wrapper function for f(t)
  use input_functions,only: f_func
  implicit none
  real(dp),intent(in) :: t
  dFdt=f_func(t)
 end function dFdt
 
 real(dp) function d2Fdt2(t)
  !!wrapper function for derivative of f(t)
  use input_functions,only: f_derivative
  implicit none
  real(dp),intent(in) :: t
  d2Fdt2=f_derivative(t)
 end function d2Fdt2
 
 real(dp) function signum(x)
  !!signum function -- equivalent to Maple's abs(1,x) function
  !!assumes Im[x]=0 and x != 0
  implicit none
  complex(dp),intent(in) :: x
  if (x%RE.lt.0._dp) then
   signum=-1._dp
  else if (x%RE.gt.0._dp) then
   signum=1._dp
  else
   signum=0._dp !!this branch is not used in practice
  end if
 end function signum
 
 real(dp) function Heaviside(x)
  !!Heavyside function
  !!To match definition used in Maple, we take Heaviside(0)=1
  !!assumes Im[x]=0
  implicit none
  complex(dp),intent(in) :: x
  if (x%RE.lt.0._dp) then
   Heaviside=0._dp
  else
   Heaviside=1._dp
  end if
 end function Heaviside
 
 real(dp) function Dirac(x)
  !!Dirac delta function -- to avoid overflow, the spike at x=0 is ignored (relevant terms go to zero as x --> 0)
  implicit none
  complex(dp),intent(in) :: x
  Dirac=0._dp
 end function Dirac
 
 real(dp) function Dirac_derivative(x)
  !!derivative of Dirac delta function -- to avoid overflow, the spike at x=0 is ignored (relevant terms go to zero as x --> 0)
  implicit none
  complex(dp),intent(in) :: x
  Dirac_derivative=0._dp
 end function Dirac_derivative
 
 complex(dp) function InverseJacobiAM(phi,k)
  !!wrapper for trigonometric form of the incomplete elliptic integral of the first kind
  !!Assumes Im{phi}=0 and Im[k]=0
  use elliptic,only: incomplete_elliptic_integrals
  implicit none
  
  !!input
  complex(dp),intent(in) :: k !!elliptic modulus
  complex(dp),intent(in) :: phi !!elliptic amplitude
  
  !!internal variables
  complex(dp) :: F,E
  real(dp)    :: m !!elliptic parameter
  
  m=k%RE**2
  call incomplete_elliptic_integrals(phi%RE,m,F,E)
  InverseJacobiAM=F
 end function InverseJacobiAM
 
 subroutine compute_JacobiSN_CN_DN(u,k,sn,cn,dn)
  !!wrapper for sn, cn, and dn functions
  !!assumes Im[k]=0
  use elliptic,only: Jacobi_elliptic_functions
  implicit none
  
  !!input
  complex(dp),intent(in) :: k !!elliptic modulus
  complex(dp),intent(in) :: u !!argument
  
  !!output 
  complex(dp),intent(out) :: sn,cn,dn
  
  !!internal variables
  real(dp) :: m !!elliptic parameter
  
  m=k%RE**2
  call Jacobi_elliptic_functions(u,m,sn,cn,dn)
 end subroutine compute_JacobiSN_CN_DN
 
 complex(dp) function JacobiSN(u,k)
  !!wrapper for sn function
  !!assumes Im[k]=0
  use elliptic,only: Jacobi_elliptic_functions
  implicit none
  
  !!input
  complex(dp),intent(in) :: k !!elliptic modulus
  complex(dp),intent(in) :: u !!argument
  
  !!internal variables
  real(dp) :: m !!elliptic parameter
  complex(dp) :: sn,cn,dn
  
  m=k%RE**2
  call Jacobi_elliptic_functions(u,m,sn,cn,dn)
  JacobiSN=sn
 end function JacobiSN
 
 complex(dp) function JacobiCN(u,k)
  !!wrapper for cn function
  !!assumes Im(k)=0
  use elliptic,only: Jacobi_elliptic_functions
  implicit none
  
  !!input
  complex(dp),intent(in) :: k !!elliptic modulus
  complex(dp),intent(in) :: u !!argument
  
  !!internal variables
  real(dp)    :: m !!elliptic parameter
  complex(dp) :: sn,cn,dn
  
  m=k%RE**2
  call Jacobi_elliptic_functions(u,m,sn,cn,dn)
  JacobiCN=cn
 end function JacobiCN
 
 complex(dp) function JacobiDN(u,k)
  !!wrapper for dn function
  !!assumes Im[k]=0
  use elliptic,only: Jacobi_elliptic_functions
  implicit none
  
  !!input
  complex(dp),intent(in) :: k !!elliptic modulus
  complex(dp),intent(in) :: u !!argument
  
  !!internal variables
  real(dp)    :: m !!elliptic parameter
  complex(dp) :: sn,cn,dn
  
  m=k%RE**2
  call Jacobi_elliptic_functions(u,m,sn,cn,dn)
  JacobiDN=dn
 end function JacobiDN
 
 complex(dp) function EllipticK(k)
  !!wrapper for complete elliptic integral of the first kind
  !!assumes Im[k]=0
  use elliptic,only: complete_elliptic_integrals
  implicit none
  
  !!input
  complex(dp),intent(in) :: k !!elliptic modulus
  
  !!internal variables
  real(dp)    :: m !!elliptic parameter
  complex(dp) :: Fc,Ec
  
  m=k%RE**2
  call complete_elliptic_integrals(m,Fc,Ec)
  EllipticK=Fc
 end function EllipticK
 
 complex(dp) function EllipticE(k)
  !!wrapper for complete elliptic integral of the second kind
  !!assumes Im[k]=0
  use elliptic,only: complete_elliptic_integrals
  implicit none
  
  !!input
  complex(dp),intent(in) :: k !!elliptic modulus
  
  !!internal variables
  real(dp)    :: m !!elliptic parameter
  complex(dp) :: Fc,Ec
  
  m=k%RE**2
  call complete_elliptic_integrals(m,Fc,Ec)
  EllipticE=Ec
 end function EllipticE
 
 subroutine compute_EllipticK_EllipticE(k,EllipticK,EllipticE)
  !!wrapper for complete elliptic integral of the first and second kinds
  !!assumes Im[k]=0
  use elliptic,only: complete_elliptic_integrals
  implicit none
  
  !!input
  complex(dp),intent(in) :: k !!elliptic modulus
  
  !!output
  complex(dp),intent(out) :: EllipticK,EllipticE
  
  !!internal variables
  real(dp)    :: m !!elliptic parameter
  complex(dp) :: Fc,Ec
  
  m=k%RE**2
  call complete_elliptic_integrals(m,Fc,Ec)
  EllipticK=Fc; EllipticE=Ec
 end subroutine compute_EllipticK_EllipticE

end module H_helper_routines
