complex*16 function sech(u)
!!hyperbolic secant function
implicit none

!!input
complex*16 :: u

!! internal variables
complex*16 :: c

c = cosh(u)
if(abs(c) > huge(0.d0)) then 
 ! For some reason 1/(inf+1i*inf) = (NaN,NaN) in Fortran, so we'll manually account for that
 sech=0.d0
else
 sech=1.d0/c
end if
end function sech

subroutine Jacobi_elliptic_functions(u,m,sn,cn,dn)
!!computes the Jacobi elliptic functions sn, cn, and dn
!!argument u can be any complex value
!!real elliptic parameter must satisfy 0<=m (m>1 is allowed)
implicit none

!!input
complex*16 :: u !!argument 
real*8 :: m !!elliptic parameter

!!output
complex*16 :: sn,cn,dn !!Jacobi elliptic functions

!!internal variables
real*8 :: mr !!reciprocal parameter
real*8 :: k  !!modulus
complex*16 :: u_temp !!temporary argument
complex*16 :: sn_temp,cn_temp,dn_temp !!temporary Jacobi elliptic function values
complex*16 :: sech

if (abs(u) < 8*epsilon(0.d0)) then
 ! sn,cn,dn are constant wrt m when u=0. This case is useful because z0 will call
 ! this subroutine with m=\infty and u\approx0 on H's boundary. In this case, u_temp=NaN
 sn=0.d0; dn=1.d0; cn=1.d0
elseif (abs(m-1.d0) < 8*epsilon(m)) then
 ! m=1 cases from fig 16.6 of Abramowitz and Stegun Handbook of Mathematical Functions
 sn=tanh(u); cn=sech(u); dn=sech(u)
elseif (m.gt.1.d0) then
 k=sqrt(m); mr=1.d0/m
 u_temp=k*u
 !write(*,*) "JEF u,m,k,u_temp,mr=",u,m,k,u_temp,mr
 call Jacobi_elliptic_functions_complex_argument_standard_parameter(u_temp,real(mr,16),sn_temp,cn_temp,dn_temp)
 sn=sn_temp/k; dn=cn_temp; cn=dn_temp
else
 call Jacobi_elliptic_functions_complex_argument_standard_parameter(u,real(m,16),sn,cn,dn)
end if
end subroutine Jacobi_elliptic_functions

subroutine Jacobi_elliptic_functions_complex_argument_standard_parameter(u,m,sn,cn,dn)
!!u can be any complex value
!!assumes 0<=m<=1
implicit none

!!input
complex*16 :: u !!argument 
real*16 :: m !!elliptic parameter

!!output
complex*16 :: sn,cn,dn !!Jacobi elliptic functions

!!internal variables
real*8 :: u_temp !!temporary argument
real*8 :: u_temp_c !!temporary argument associated with complimentary parameter
real*8 :: sn_temp,cn_temp,dn_temp !!temporary Jacobi elliptic
real*8 :: sn_temp_c,cn_temp_c,dn_temp_c !!temporary Jacobi elliptic using complimentary parameters
real*16 :: mc !!compliment of parameter
real*8 :: delta !!divisor
real*8 :: sn_temp_dn_temp_c,cn_temp_cn_temp_c,dn_temp_cn_temp_c !!products of temporary variables
real*8 :: k
real*8 :: m_dp !!double precision value of m

if (u%IM.ne.0.d0) then !!complex argument
 m_dp=real(m,8)
 k=sqrt(m_dp)
 mc=1.q0-m
 !mc=(1.d0+k)*(1.d0-k)
 u_temp_c=u%IM
 !write(*,*) "JEFCargSParam: u,m=",u,m
 !write(*,*) "JEFCargSParam: u_temp_c,mc=",u_temp_c,mc
 call Jacobi_elliptic_functions_standard_input_range(u_temp_c,mc,sn_temp_c,cn_temp_c,dn_temp_c)
 !write(*,*) "JEFCargSParam: sn_temp_c,cn_temp_c,dn_temp_c=",sn_temp_c,cn_temp_c,dn_temp_c
 !write(*,*) ""
 u_temp=u%RE
 !write(*,*) "JEFCargSParam: u_temp,m=",u_temp,m
 call Jacobi_elliptic_functions_standard_input_range(u_temp,m,sn_temp,cn_temp,dn_temp)
 !write(*,*) "JEFCargSParam: sn_temp,cn_temp,dn_temp=",sn_temp,cn_temp,dn_temp
 delta=cn_temp_c**2.d0+m_dp*(sn_temp*sn_temp_c)**2.d0 !!reduces truncation error (this form avoids subtraction)

 sn_temp_dn_temp_c=sn_temp*dn_temp_c; cn_temp_cn_temp_c=cn_temp*cn_temp_c; dn_temp_cn_temp_c=dn_temp*cn_temp_c

 sn=complex(sn_temp_dn_temp_c,sn_temp_c*cn_temp_cn_temp_c*dn_temp)/delta
 cn=complex(cn_temp_cn_temp_c,-sn_temp_dn_temp_c*dn_temp*sn_temp_c)/delta
 dn=complex(dn_temp_cn_temp_c*dn_temp_c,-m_dp*sn_temp*cn_temp*sn_temp_c)/delta

else !!real argument
 u_temp=u%RE
 call Jacobi_elliptic_functions_standard_input_range(u_temp,m,sn_temp,cn_temp,dn_temp)
 sn=complex(sn_temp,0.d0); cn=complex(cn_temp,0.d0); dn=complex(dn_temp,0.d0)
end if
end subroutine Jacobi_elliptic_functions_complex_argument_standard_parameter

subroutine Jacobi_elliptic_functions_standard_input_range(u,m,sn,cn,dn)
!!computes the Jacobi elliptic functions sn, cn, and dn
!!argument u can be any real value
!!elliptic parameter must be within 0<=m<=1  **** m is quadruple precision ****
implicit none

!!input
real*8 :: u !!argument and elliptic parameter
real*16 :: m

!!output
real*8 :: sn,cn,dn !!Jacobi elliptic functions

!!internal variables
real*16 :: mc !!compliment of elliptic parameter
!real*8 :: k !!elliptic modulus

!k=sqrt(m)
!mc=(1.d0+k)*(1.d0-k)
mc=1.q0-m
call gscd(u,mc,sn,cn,dn)
end subroutine Jacobi_elliptic_functions_standard_input_range

subroutine complete_elliptic_integrals(m,Fc,Ec)
!!!!SJT: computes the complete elliptic integrals of first and second kind given parameter m
!!!!SJT: assumes 0<=m
!!!!SJT: allows m>1 by applying the reciprocal-modulus transformation
implicit none

!!input
real*8 :: m !!elliptic characteristic and parameter

!!output
complex*16 :: Fc,Ec !!complete elliptic integrals of the first and second kinds

!!internal variables
real*8 :: mc !!compliment of elliptic parameter
real*8 :: mr,mrc !!reciprocal parameter and its compliment 
real*8 :: k !!elliptic modulus and reciprocal
real*8 :: n !!characteristic
real*8 :: Fc_temp,Ec_temp,Pc_temp !!temporary complete elliptic integral values from standard input range
real*8 :: Fc_r,Ec_r,Pc_r,Fc_rc,Ec_rc,Pc_rc !!complete elliptic integral values based on reciprocal parameter and its compliment
real*16 :: m_qp,mc_qp,mr_qp,mrc_qp !!quad precision variables

n=0.d0 !!arbitrary characteristic value -- integrals of the third kind are not used

if (m.gt.1.d0) then
 ! Note that when m is large (eg m>=10^30) then mr_qp is approx 0 and mrc_qp is approx 1
 ! which can be problematic
 k=sqrt(m)
 m_qp=real(m,16); mr_qp=real(1.d0/m,16); mrc_qp=1.q0-mr_qp; mc_qp=1.q0-m_qp; mc=real(mc_qp,8)
 call complete_elliptic_integrals_standard_input_range(n,mr_qp,Fc_r,Ec_r,Pc_r)
 call complete_elliptic_integrals_standard_input_range(n,mrc_qp,Fc_rc,Ec_rc,Pc_rc)
 Fc=complex(Fc_r,-Fc_rc)/k; Ec=complex(m*Ec_r+mc*Fc_r,-Fc_rc+m*Ec_rc)/k !!OG -- cancellation error possible for real part of Ec
else
 call complete_elliptic_integrals_standard_input_range(n,real(m,16),Fc_temp,Ec_temp,Pc_temp)
 Fc=complex(Fc_temp,0.d0); Ec=complex(Ec_temp,0.d0)
end if
end subroutine complete_elliptic_integrals

subroutine complete_elliptic_integrals_standard_input_range(n,m,Fc,Ec,Pc)
!!!!SJT: computes the complete elliptic integrals of kinds first to third given parameter m and characteristic n
!!!!SJT: assumes 0<=m<=1 and 0<=n<=1
implicit none

!!input
real*8 :: n !!elliptic characteristic and parameter
real*16 :: m

!!output
real*8 :: Fc,Ec,Pc !!complete elliptic integrals of the first kind, second kind, and third kind

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0
real*8 :: zero,piio2 !!constants
real*8 :: bc,dc,jc !!associate complete elliptic integrals
real*16 :: mc !!complimentary parameter
!real*8 :: k !!elliptic modulus

zero=0.d0
piio2=pii/2.d0

!mc=1.d0-m !!complimentary parameter
mc=1.q0-m
!k=sqrt(m); mc=(1.d0+k)*(1.d0-k)

! Fc->\infty EXTREMELY SUDDENLY as m->1. For example, when m=1-1e-33 then Fc\approx 40.
! Thankfully, elbdj2 only returns NaNs when m is literally exactly 1. It'd likely converge
! for more extreme values of m if the type of m could represent those values
if (m.eq.1.d0) then ! strict equality is intentional here
 ! when m=1 the integrals defining Fc,Ec,Pc are straightforward to compute with a trig sub
 Fc = 50.d0; Ec = 1.d0
 !Fc = huge(0.d0); Ec = 1.d0
 if (abs(n) < 4*epsilon(n)) then
  jc = 1/3.d0
 else if (n.ge.1) then
  jc = huge(0.d0) ! TODO this might be too naive... maybe there's a well defined complex output?
 else if (n > 0) then
  jc = atanh(sqrt(n)*sin(piio2))/(n**1.5) - sin(piio2)/n
 else ! n < 0
  jc = -sin(piio2)/n - atan(sqrt(-n)*sin(piio2))/((-n)**1.5)
 end if
 Pc=Fc+n*jc
else
 call elbdj2(piio2,zero,n,mc,bc,dc,jc) !!returns the associate complete elliptic integrals
 Fc=bc+dc; Ec=bc+real(mc,8)*dc; Pc=Fc+n*jc !!build standard elliptic integrals from associate elliptic integrals
end if
end subroutine complete_elliptic_integrals_standard_input_range

subroutine incomplete_elliptic_integrals_standard_input_range(phi,n,m,F,E,P)
!!!!SJT: computes the incomplete elliptic integrals of kinds first to third given parameter m and characteristic n
!!assumes 0<=phi<=pi/2, 0<=m<=1, 0<=n<=1
implicit none

!!input
real*8 :: phi,n !!elliptic amplitude, characteristic, and parameter
real*16 :: m

!!output
real*8 :: F,E,P !!incomplete elliptic integrals of the first kind, second kind, and third kind

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0
real*8 :: piio2      !!constants
real*8 :: phic    !!complimentary variables
real*16 :: mc
!real*8 :: k          !!elliptic modulus
real*8 :: b,d,j      !!associate incomplete elliptic integrals

piio2=pii/2.d0

phic=piio2-phi;
!k=sqrt(m); mc=(1.d0+k)*(1.d0-k)
!mc=1.d0-m
mc=1.q0-m
if (abs(m-1.d0).eq.epsilon(m)) then !! parameter m is approx 1
 F = log(1.d0/cos(phi) + tan(phi)); E = sign(sin(phi), cos(phi))
 if (abs(n) < 8*epsilon(n)) then
  j = sin(phi)**3/3
 else if (n>0) then
  j = atanh(sqrt(n)*sin(phi))/(n**1.5) - sin(phi)/n
 else ! n < 0
  j = -sin(phi)/n - atan(sqrt(-n)*sin(phi))/((-n)**1.5)
 end if
 P=F+n*j
else
 call elbdj2(phi,phic,n,mc,b,d,j) !!associate incomplete elliptic integrals
 F=b+d; E=b+real(mc,8)*d; P=F+n*j !!build standard elliptic integrals from associate elliptic integrals
end if
end subroutine incomplete_elliptic_integrals_standard_input_range

subroutine incomplete_elliptic_integrals_standard_amp_large_parameter(phi,m,F,E)
!!!!SJT: reduce parameter using reciprocal-modulus transformation (if needed)
!!!!SJT: computes the incomplete elliptic integrals of first and second kind given amplitude phi and parameter m
!!assumes 0<=phi<=pi/2 and 0<=m (m>1 is allowed)
implicit none

!!input
real*8 :: phi,m !!elliptic amplitude and parameter

!!output
complex*16 :: F,E !!incomplete elliptic integrals of the first and second kind

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0
real*16, parameter :: pii_qp=3.14159265358979323846264338327950288q0
real*8 :: piio2 !!constants
real*8 :: phic,mc    !!complimentary variables
real*8 :: phi_temp,phic_temp,b_temp,d_temp,j_temp !!temporary variables
real*8 :: u !!argument of arcsine
real*8 :: k !!elliptic modulus
real*8 :: kr,mr !!reciprocal modulus and parameter
real*8 :: mrc !!compliment of the reciprocal parameter
complex*16 :: Fc,Ec !!complete elliptic integral values computed for m>1
real*8 :: Fc_temp,Ec_temp,Pc_temp !!temporary complete elliptic integral values
real*8 :: F_temp,E_temp,P_temp !!temporary incomplete elliptic integral values
real*8 :: amp,ampc !!temporary elliptic amplitude variable and compliment
real*8 :: sin_amp !!sine amplitude
real*8 :: n !!characteristic
real*16 :: mr_qp,mrc_qp,m_qp,mc_qp
real*16 :: u_qp,phi_temp_qp,amp_qp,sin_amp_qp

piio2=pii/2.d0
n=0.d0 !!arbitrary characteristic -- integrals of the third kind are not used

if (phi.eq.piio2) then
 call complete_elliptic_integrals(m,F,E)
else !!0<=phi<pi/2
 if (m.gt.1.d0) then !!large parameter m>1
  k=sqrt(m)
  !mc=1.d0-m; mr=1.d0/m !!OG
  mc=real(1.q0-real(m,16),8); mr_qp=1.q0/real(m,16)
  !u=k*sin(phi) !!argument of arcsine
  u_qp=real(k,16)*sin(real(phi,16))
  if (u_qp.le.1.q0) then
   !amp=asin(u)
   amp=real(asin(u_qp),8)
   call incomplete_elliptic_integrals_standard_input_range(amp,n,mr_qp,F_temp,E_temp,P_temp)
   F=complex(F_temp/k,0.d0); E=complex(k*E_temp+mc*F%RE,0.d0)
  else
   !mrc=1.d0-mr; !!OG
   mrc_qp=1.q0-mr_qp; mrc=real(mrc_qp,8)
   !sin_amp=(sqrt(u**2.d0-1.d0))/(u*sqrt(mrc))
   sin_amp_qp=(sqrt(u_qp**2-1.q0))/(u_qp*sqrt(mrc_qp))
   !amp=asin(min(sin_amp,1.d0)) !!it is possible for sin_amp to be slightly greater than unity due to round off error
   amp_qp=asin(min(sin_amp_qp,1.q0)); amp=real(amp_qp,8) !!it may be possible for sin_amp to be slightly greater than unity due to round off error
   call complete_elliptic_integrals_standard_input_range(n,mr_qp,Fc_temp,Ec_temp,Pc_temp)
   call incomplete_elliptic_integrals_standard_input_range(amp,n,mrc_qp,F_temp,E_temp,P_temp)
   F=complex(Fc_temp,-F_temp)/k
   E=complex(m*Ec_temp+mc*Fc_temp,&
     &-F_temp+m*E_temp+real(((1.q0-m)*sin_amp_qp*cos(amp_qp))/sqrt(1.q0-mrc_qp*sin_amp_qp**2),8))/k
  end if
 else !!standard input ranges
  m_qp=real(m,16)
  call incomplete_elliptic_integrals_standard_input_range(phi,n,m_qp,F_temp,E_temp,P_temp)
  F=complex(F_temp,0.d0); E=complex(E_temp,0.d0)
 end if
end if

end subroutine incomplete_elliptic_integrals_standard_amp_large_parameter

subroutine incomplete_elliptic_integrals(phi,m,F,E)
!!!!SJT: computes the incomplete elliptic integrals of first and second kind given amplitude phi and parameter m
!!assumes phi is real
!!assumes 0<=m (m>1 is allowed but m must be real)
implicit none

!!input
real*8 :: phi,m !!elliptic amplitude and parameter

!!output
complex*16 :: F,E !!incomplete elliptic integrals of the first and second kind

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0
real*8 :: piio2 !!constants
complex*16 :: Fc,Ec !!complete elliptic integral values computed for m>1
real*8 :: phi_std
real*8 :: N0
integer*4 :: N,i_sign

piio2=pii/2.d0

if ((0.d0.le.phi).and.(phi.le.piio2)) then !!standard amplitude range
 call incomplete_elliptic_integrals_standard_amp_large_parameter(phi,m,F,E)
else !!amplitude outside of standard range
 call complete_elliptic_integrals(m,Fc,Ec)
 N0=phi/pii;
 N=nint(N0,4) !!integer multiplier
 
 phi_std=phi-N*pii !!find appropriate phi in standard range
 i_sign = sign(1.d0,phi_std)
 phi_std=abs(phi_std)

 call incomplete_elliptic_integrals_standard_amp_large_parameter(phi_std,m,F,E)
 F=2*N*Fc+i_sign*F
 E=2*N*Ec+i_sign*E
end if

end subroutine incomplete_elliptic_integrals

