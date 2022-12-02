complex*16 function JacobiZeta(u,k)
!!Jacobi Zeta function
!!assumes u belongs to [F(0|m),F(pi|m)]
!!assumes imaginary part of k is zero
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument (equivalent to F(phi|m))

!!internal variables
real*8 :: m !!elliptic parameter
real*8 :: phi !!Jacobi amplitude
complex*16 :: Fc,Ec,F,E !!elliptic integral values
complex*16 :: sn,cn,dn !!Jacobi elliptic function values

!!!original
!m=k%RE**2
!call complete_elliptic_integrals(m,Fc,Ec)
!call Jacobi_elliptic_functions(u,m,sn,cn,dn)
!
!phi=acos(cn%RE) !!due to the presumed range of u
!call incomplete_elliptic_integrals(phi,m,F,E) !!assuming imaginary part of phi is negligible
!
!JacobiZeta=E-Ec*u/Fc
!!end original

call compute_JacobiZeta(u,k,JacobiZeta)
end function JacobiZeta

subroutine compute_JacobiZeta(u,k,JacobiZeta)
!!Jacobi Zeta function
!!assumes u belongs to [F(0|m),F(pi|m)]
!!assumes imaginary part of k is zero
!!assumes the Jacobi amplitude (phi) is real
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument (equivalent to F(phi|m))

!!output
complex*16 :: JacobiZeta

!!internal variables
real*8, parameter :: tol=1.d-7 !!base tolerance for integer tests
real*8 :: m !!elliptic parameter
real*8 :: phi !!Jacobi amplitude
complex*16 :: Fc,Ec,F,E !!elliptic integral values
complex*16 :: twoFc !!two times Fc
complex*16 :: sn,cn,dn !!Jacobi elliptic function values
real*8 :: Nm_IM,Np_IM,Nm_RE,Np_RE !!periodicity multiplier
real*8 :: fraction_Nm_RE,fraction_Np_RE,fraction_Nm_IM,fraction_Np_IM
logical :: IM_periodic,RE_periodic
logical :: int_p_RE,int_m_RE,int_p_IM,int_m_IM

m=k%RE**2
call complete_elliptic_integrals(m,Fc,Ec)
call Jacobi_elliptic_functions(u,m,sn,cn,dn)

if (abs(cn%RE).gt.1.d0) then
 write(*,*) "Error in compute_JacobiZeta -- abs(cn%RE) is greater than unity"
 write(*,*) "cn=",cn
 stop
elseif (abs(cn%IM).gt.1.d-10) then
 write(*,*) "Error in compute_JacobiZeta -- cn%IM is not zero"
 write(*,*) "u,k,m=",u,k,m
 write(*,*) "cn=",cn
 stop
end if
phi=acos(cn%RE) !!produces phi belonging to [0,pi] which works with incomplete_elliptic_integrals routine
call incomplete_elliptic_integrals(phi,m,F,E) !!assuming imaginary part of phi is negligible
twoFc=2.d0*Fc !!store this for frequent use

!!Re[u] not in the presumed range? 
RE_periodic=.not.((((0.d0.le.u%RE).and.(u%RE.le.(twoFc%RE))).or.(((twoFc%RE).lt.u%RE).and.(u%RE.le.0.d0))))
if (RE_periodic.eqv..true.) then
 Nm_RE=(u%RE-F%RE)/(twoFc%RE)
 Np_RE=(u%RE+F%RE)/(twoFc%RE)
 fraction_Nm_RE=abs(mod(Nm_RE,1.d0)) !fractional parts of N values -- should be very close to zero for integer values of N
 fraction_Np_RE=abs(mod(Np_RE,1.d0))
 int_p_RE=((fraction_Np_RE.le.abs(tol*Np_RE)).or.((1.d0-fraction_Np_RE).le.abs(tol*Np_RE)).or.(abs(Np_RE).lt.tol))
 int_m_RE=((fraction_Nm_RE.le.abs(tol*Nm_RE)).or.((1.d0-fraction_Nm_RE).le.abs(tol*Nm_RE)).or.(abs(Nm_RE).lt.tol))
 if (int_p_RE.eqv..true.) then
  E%RE=2.d0*Ec%RE*Np_RE-E%RE !!use periodicity of elliptic integrals
 elseif (int_m_RE.eqv..true.) then
  E%RE=2.d0*Ec%RE*Nm_RE+E%RE
 else
  write(*,*) "Error in compute_JacobiZeta -- integer multipliers for real parts do not meet tolerance."
  write(*,*) "  Either Nm_RE or Np_RE should be near an integer value."
  write(*,*) "  Nm_RE,Np_RE=",Nm_RE,Np_RE
  write(*,*) "  May need to increase tolerance (tol parameter) in compute_JacobiZeta routine (H_helper_functions.f90)."
  stop
 end if
end if

!!Im[u] not in the presumed range?
IM_periodic=.not.((((0.d0.le.u%IM).and.(u%IM.le.(twoFc%IM))).or.(((twoFc%IM).lt.u%IM).and.(u%IM.le.0.d0))))
if (IM_periodic.eqv..true.) then
 Nm_IM=(u%IM-F%IM)/(twoFc%IM)
 Np_IM=(u%IM+F%IM)/(twoFc%IM)
 fraction_Nm_IM=abs(mod(Nm_IM,1.d0)) !!fractional parts of N values -- should be very close to zero for integer values of N
 fraction_Np_IM=abs(mod(Np_IM,1.d0))
 int_p_IM=((fraction_Np_IM.le.abs(tol*Np_IM)).or.((1.d0-fraction_Np_IM).le.abs(tol*Np_IM)).or.(abs(Np_IM).lt.tol))
 int_m_IM=((fraction_Nm_IM.le.abs(tol*Nm_IM)).or.((1.d0-fraction_Nm_IM).le.abs(tol*Nm_IM)).or.(abs(Nm_IM).lt.tol))
 if (int_p_IM.eqv..true.) then
  E%IM=2.d0*Ec%IM*Np_IM-E%IM
 elseif (int_m_IM.eqv..true.) then
  E%IM=2.d0*Ec%IM*Nm_IM+E%IM
 else
  write(*,*) "Error in compute_JacobiZeta -- integer multipliers for imaginary parts do not meet tolerance."
  write(*,*) "  Either Nm_IM or Np_IM should be near an integer value."
  write(*,*) "  Nm_IM,Np_IM=",Nm_IM,Np_IM
  write(*,*) "  May need to increase tolerance (tol parameter) in compute_JacobiZeta routine (H_helper_functions.f90)."
  stop
 end if
end if

JacobiZeta=E-Ec*u/Fc
end subroutine compute_JacobiZeta

real*8 function F(t)
!!wrapper function for integral of f(t)
implicit none
real*8 :: t
real*8 :: f_integral
F=f_integral(t)
end function F

real*8 function dFdt(t)
!!wrapper function for f(t)
implicit none
real*8 :: t
real*8 :: f_func
dFdt=f_func(t)
end function dFdt

real*8 function d2Fdt2(t)
!!wrapper function for derivative of f(t)
implicit none
real*8 :: t
real*8 :: f_derivative
d2Fdt2=f_derivative(t)
end function d2Fdt2

real*8 function d3Fdt3(t)
!!wrapper function for derivative of f(t)
implicit none
real*8 :: t
real*8 :: f_2nd_derivative
d3Fdt3=f_2nd_derivative(t)
end function d3Fdt3

real*8 function Heaviside(x)
!!Heavyside function
!!To match definition used in Maple, we take Heaviside(0)=1
!!assumes Im[x]=0
implicit none
complex*16 :: x
!!if (x.le.0.d0) then
if (x%RE.lt.0.d0) then
 Heaviside=0.d0
else
 Heaviside=1.d0
end if
end function Heaviside

real*8 function Dirac(x)
!!Dirac delta function -- to avoid overflow the spike at x=0 is ignored (relevant terms go to zero as x --> 0)
implicit none
complex*16 :: x
Dirac=0.d0
end function Dirac

real*8 function Dirac_derivative(x)
!!derivative of Dirac delta function -- to avoid overflow the spike at x=0 is ignored (relevant terms go to zero as x --> 0)
implicit none
complex*16 :: x
Dirac_derivative=0.d0
end function Dirac_derivative

complex*16 function InverseJacobiAM(phi,k)
!!wrapper for trigonometric form of the incomplete elliptic integral of the first kind
!!Assumes Im{phi}=0 and Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: phi !!elliptic amplitude

!!internal variables
complex*16 :: F,E
real*8 :: m !!elliptic parameter

m=k%RE**2
call incomplete_elliptic_integrals(phi%RE,m,F,E)
InverseJacobiAM=F
end function InverseJacobiAM

subroutine compute_JacobiSN_CN_DN(u,k,sn,cn,dn)
!!wrapper for sn, cn, and dn functions
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument

!!output 
complex*16 :: sn,cn,dn

!!internal variables
real*8 :: m !!elliptic parameter

m=k%RE**2
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
end subroutine compute_JacobiSN_CN_DN

complex*16 function JacobiSN(u,k)
!!wrapper for sn function
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: sn,cn,dn

m=k%RE**2
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
JacobiSN=sn
end function JacobiSN

complex*16 function JacobiCN(u,k)
!!wrapper for cn function
!!assumes Im(k)=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: sn,cn,dn

m=k%RE**2
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
JacobiCN=cn
end function JacobiCN

complex*16 function JacobiDN(u,k)
!!wrapper for dn function
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: sn,cn,dn

m=k%RE**2
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
JacobiDN=dn
end function JacobiDN

complex*16 function EllipticK(k)
!!wrapper for complete elliptic integral of the first kind
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: Fc,Ec

m=k%RE**2
call complete_elliptic_integrals(m,Fc,Ec)
EllipticK=Fc
end function EllipticK

complex*16 function EllipticE(k)
!!wrapper for complete elliptic integral of the second kind
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: Fc,Ec

m=k%RE**2
call complete_elliptic_integrals(m,Fc,Ec)
EllipticE=Ec
end function EllipticE

subroutine compute_EllipticK_EllipticE(k,EllipticK,EllipticE)
!!wrapper for complete elliptic integral of the first and second kinds
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus

!!output
complex*16 :: EllipticK,EllipticE

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: Fc,Ec

m=k%RE**2
call complete_elliptic_integrals(m,Fc,Ec)
EllipticK=Fc; EllipticE=Ec
end subroutine compute_EllipticK_EllipticE
