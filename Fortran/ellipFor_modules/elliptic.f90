module elliptic
 ! module containing main ellipFor subroutines for the calculation of Legendre elliptic integrals 
 ! and Jacobi elliptic functions with generalized input ranges
 use kind_parameters
 implicit none
 private
 public :: Jacobi_elliptic_functions     ! compute the primary Jacobi elliptic functions sn, cn, and dn
 public :: complete_elliptic_integrals   ! compute complete Legendre elliptic intergrals of the first and second kinds
 public :: incomplete_elliptic_integrals ! compute incomplete Legendre elliptic intergrals of the first and second kinds

 ! parameters
 real(dp), parameter :: pii   =3.1415926535897932_dp
 real(qp), parameter :: pii_qp=3.14159265358979323846264338327950288_qp
contains
 
 elemental subroutine Jacobi_elliptic_functions(u,m,sn,cn,dn)
 !!computes the Jacobi elliptic functions sn, cn, and dn
 !!argument u can be any complex value
 !!real elliptic parameter must satisfy m>=0 (otherwise NaN is returned)
 use,intrinsic :: ieee_arithmetic,only: ieee_value,ieee_quiet_nan
 
 !!input
 complex(dp),intent(in) :: u !!argument 
 real(dp),intent(in) :: m !!elliptic parameter
 
 !!output
 complex(dp),intent(out) :: sn,cn,dn !!Jacobi elliptic functions
 
 !!internal variables
 real(qp)    :: mr                      !!reciprocal parameter
 real(dp)    :: k                       !!modulus
 complex(dp) :: u_temp                  !!temporary argument
 complex(dp) :: sn_temp,cn_temp,dn_temp !!temporary Jacobi elliptic function values
 real(dp)    :: nan                     !!NaN value
 real(qp)    :: m_qp                    !!quad precision m  

 if (m.gt.1._dp) then
  m_qp=real(m,qp)
  mr=1._qp/m_qp
  k=real(sqrt(m_qp),dp)
  !k=sqrt(m); mr=1.d0/m
  u_temp=k*u
  call Jacobi_elliptic_functions_complex_argument_standard_parameter(u_temp,mr,sn_temp,cn_temp,dn_temp)
  sn=sn_temp/k; dn=cn_temp; cn=dn_temp
 else if (m.gt.0._dp) then
  m_qp=real(m,qp)
  call Jacobi_elliptic_functions_complex_argument_standard_parameter(u,m_qp,sn,cn,dn)
 else if (m.eq.0._dp) then ! special case of m=0
  sn=sin(u); cn=cos(u); dn=1._dp
 else ! m < 0 not currently supported -- returning NaN values
  nan=ieee_value(0._dp,ieee_quiet_nan)
  sn=cmplx(nan,nan,dp); cn=cmplx(nan,nan,dp); dn=cmplx(nan,nan,dp) 
 end if
 end subroutine Jacobi_elliptic_functions
 
 elemental subroutine Jacobi_elliptic_functions_complex_argument_standard_parameter(u,m,sn,cn,dn)
 !!u can be any complex value
 !!assumes 0<=m<=1
 
 !!input
 complex(dp),intent(in) :: u !!argument 
 real(qp),intent(in) :: m !!elliptic parameter
 
 !!output
 complex(dp),intent(out) :: sn,cn,dn !!Jacobi elliptic functions
 
 !!internal variables
 real(dp) :: u_temp !!temporary argument
 real(dp) :: u_temp_c !!temporary argument corresponding to complimentary parameter
 real(dp) :: sn_temp,cn_temp,dn_temp !!temporary Jacobi elliptic
 real(dp) :: sn_temp_c,cn_temp_c,dn_temp_c !!temporary Jacobi elliptic using complimentary parameters
 real(qp) :: mc !!compliment of parameter
 real(dp) :: delta !!divisor
 real(dp) :: sn_temp_dn_temp_c,cn_temp_cn_temp_c,dn_temp_cn_temp_c !!products of temporary variables
 real(dp) :: m_dp !!double precision value of m
 
 if (u%IM.ne.0._dp) then !!complex argument
  m_dp=real(m,dp)
  mc=1.0_qp-m
  u_temp_c=u%IM
  call Jacobi_elliptic_functions_standard_input_range(u_temp_c,mc,sn_temp_c,cn_temp_c,dn_temp_c)
  u_temp=u%RE
  call Jacobi_elliptic_functions_standard_input_range(u_temp,m,sn_temp,cn_temp,dn_temp)
  delta=cn_temp_c**2_isp+m_dp*(sn_temp*sn_temp_c)**2_isp !!reduces truncation error (this form avoids subtraction)
 
  sn_temp_dn_temp_c=sn_temp*dn_temp_c; cn_temp_cn_temp_c=cn_temp*cn_temp_c; dn_temp_cn_temp_c=dn_temp*cn_temp_c
 
  sn=cmplx(sn_temp_dn_temp_c,sn_temp_c*cn_temp_cn_temp_c*dn_temp,dp)/delta
  cn=cmplx(cn_temp_cn_temp_c,-sn_temp_dn_temp_c*dn_temp*sn_temp_c,dp)/delta
  dn=cmplx(dn_temp_cn_temp_c*dn_temp_c,-m_dp*sn_temp*cn_temp*sn_temp_c,dp)/delta
 
 else !!real argument
  u_temp=u%RE
  call Jacobi_elliptic_functions_standard_input_range(u_temp,m,sn_temp,cn_temp,dn_temp)
  sn=cmplx(sn_temp,0._dp,dp); cn=cmplx(cn_temp,0._dp,dp); dn=cmplx(dn_temp,0._dp,dp)
 end if
 end subroutine Jacobi_elliptic_functions_complex_argument_standard_parameter
 
 elemental subroutine Jacobi_elliptic_functions_standard_input_range(u,m,sn,cn,dn)
 !!computes the Jacobi elliptic functions sn, cn, and dn
 !!argument u can be any real value
 !!elliptic parameter must be within 0<=m<=1  **** m is quadruple precision ****
 use xgscd_routines, only: gscd 
 
 !!input
 real(dp),intent(in) :: u !!argument and elliptic parameter
 real(qp),intent(in) :: m
 
 !!output
 real(dp),intent(out) :: sn,cn,dn !!Jacobi elliptic functions
 
 !!internal variables
 real(qp) :: mc !!compliment of elliptic parameter
 
 mc=1.0_qp-m
 call gscd(u,mc,sn,cn,dn)
 end subroutine Jacobi_elliptic_functions_standard_input_range
 
 elemental subroutine complete_elliptic_integrals(m,Fc,Ec)
 !!!!SJT: computes the complete elliptic integrals of first and second kind given parameter m
 !!!!SJT: assumes m>=0 (returns NaN otherwise)
 !!!!SJT: allows m>1 by applying the reciprocal-modulus transformation for complete elliptic integrals
 use,intrinsic :: ieee_arithmetic,only: ieee_value,ieee_quiet_nan,ieee_positive_inf
 
 !!input
 real(dp),intent(in) :: m         !!elliptic characteristic and parameter
 
 !!output
 complex(dp),intent(out) :: Fc,Ec !!complete elliptic integrals of the first and  second kinds
 
 !!internal variables
 real(dp) :: mc                               !!compliment of elliptic parameter
 real(dp) :: k                                !!elliptic modulus and reciprocal
 real(dp) :: n                                !!characteristic
 real(dp) :: Fc_temp,Ec_temp,Pc_temp          !!temporary complete elliptic integral values from standard input range
 real(dp) :: Fc_r,Ec_r,Pc_r,Fc_rc,Ec_rc,Pc_rc !!complete elliptic integral values based on reciprocal parameter and compliment
 real(dp) :: bc_r,dc_r,jc_r                   !!associated complete elliptic integrals
 real(dp) :: bc_rc,dc_rc,jc_rc                !!associated complete elliptic integrals
 real(dp) :: bc_temp,dc_temp,jc_temp          !!associated complete elliptic integrals
 real(qp) :: m_qp,mc_qp,mr_qp,mrc_qp          !!quad precision variables
 real(dp) :: nan                              !!NaN value
 real(dp) :: pos_inf                          !!positive infinity
 
 n=0._dp !!arbitrary characteristic value -- integrals of the third kind are not used
 
 if (m.gt.1._dp) then
  k=sqrt(m)
  m_qp=real(m,qp); mr_qp=real(1._dp/m,qp); mrc_qp=1.0_qp-mr_qp; mc_qp=1.0_qp-m_qp; mc=real(mc_qp,dp)
  call complete_elliptic_integrals_standard_input_range(n,mr_qp,Fc_r,Ec_r,Pc_r,bc_r,dc_r,jc_r)
  call complete_elliptic_integrals_standard_input_range(n,mrc_qp,Fc_rc,Ec_rc,Pc_rc,bc_rc,dc_rc,jc_rc)
  !!OG -- cancellation error possible for Re[Ec]
  !Fc=cmplx(Fc_r,-Fc_rc,dp)/k; Ec=cmplx(m*Ec_r+mc*Fc_r,-Fc_rc+m*Ec_rc,dp)/k
 
  !!reduce cancellation error by factor of 10 (but cancellation in Ec_r-Fc_r still possible)
  !!note: mc=1-m identity used
  !Fc=cmplx(Fc_r,-Fc_rc,dp)/k; Ec=cmplx(m*(Ec_r-Fc_r)+Fc_r,-Fc_rc+m*Ec_rc,dp)/k

  !!eliminate cancellation error by using identity with associated integral dc_r: Ec_r-Fc_r=-m*dc_r
  !Fc=cmplx(Fc_r,-Fc_rc,dp)/k; Ec=cmplx(m*(-real(mr_qp,dp)*dc_r)+Fc_r,-Fc_rc+m*Ec_rc,dp)/k

  !!simplify using m*real(mr_qp,dp)=1._dp (no cancellation error)
  !Fc=cmplx(Fc_r,-Fc_rc,dp)/k; Ec=cmplx(-dc_r+Fc_r,-Fc_rc+m*Ec_rc,dp)/k

  !!simplify using -dc_r+Fc_r=bc_r (eliminate differences completely in Re[Ec])
  !!note: Re[Ec] may agree better with Mathematica than SageMath for very large m  
  !Fc=cmplx(Fc_r,-Fc_rc,dp)/k; Ec=cmplx(bc_r,-Fc_rc+m*Ec_rc,dp)/k

  !!express Im[Ec] in terms of bc_rc to reduce error near m=1
  !!note: used Fc_rc=bc_rc+dc_rc, Ec_rc=bc_rc+(1-mrc)*dc_rc, 1-mrc=1/m
  Fc=cmplx(Fc_r,-Fc_rc,dp)/k; Ec=cmplx(bc_r,-mc*bc_rc,dp)/k
 else if (m.eq.1._dp) then ! special values
  pos_inf = ieee_value(0._dp,ieee_positive_inf) ! +Inf
  Fc=cmplx(pos_inf,pos_inf,dp)                  ! complex infinity: no ieee value so using (+Inf,+Inf)
  Ec=cmplx(1._dp,0._dp,dp)                      ! unity
 else if (m.ge.0._dp) then
  call complete_elliptic_integrals_standard_input_range(n,real(m,qp),Fc_temp,Ec_temp,Pc_temp,bc_temp,dc_temp,jc_temp)
  Fc=cmplx(Fc_temp,0._dp,dp); Ec=cmplx(Ec_temp,0._dp,dp)
 else ! m < 0 not currently supported -- returning NaN values
  nan=ieee_value(0._dp,ieee_quiet_nan)
  Fc=cmplx(nan,nan,dp); Ec=cmplx(nan,nan,dp)
 end if
 end subroutine complete_elliptic_integrals
 
 elemental subroutine complete_elliptic_integrals_standard_input_range(n,m,Fc,Ec,Pc,bc,dc,jc)
 !!!!SJT: computes the complete elliptic integrals of kinds first to third given parameter m and characteristic n
 !!!!SJT: assumes 0<=m<=1 and 0<=n<=1
 use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_positive_inf
 use xelbdj2_all_routines, only: elbdj2 

 !!input
 real(dp),intent(in) :: n !!elliptic characteristic and parameter
 real(qp),intent(in) :: m
 
 !!output
 real(dp),intent(out) :: Fc,Ec,Pc !!complete elliptic integrals of the first kind, second kind, and third kind
 real(dp),intent(out) :: bc,dc,jc !!associated complete elliptic integrals
 
 !!internal variables
 real(dp) :: zero,piio2 !!constants
 real(qp) :: mc         !!complimentary parameter

 if (m.eq.1._dp) then
  Fc=ieee_value(1._dp,ieee_positive_inf)
  Ec=1._dp
  Pc=0._dp !! arbitrarily set to zero for now -- update in future ellipFor versions
 else 
  zero=0._dp
  piio2=pii/2._dp
  
  mc=1.0_qp-m !!complimentary parameter (quad precision)
  call elbdj2(piio2,zero,n,mc,bc,dc,jc) !!returns the associated complete elliptic integrals
  Fc=bc+dc; Ec=bc+real(mc,dp)*dc; Pc=Fc+n*jc !!OG: build standard elliptic integrals from associated elliptic integrals
  !Fc=bc+dc; Ec=Fc-real(m,dp)*dc; Pc=Fc+n*jc !!slightly worse for Ec (but useful identity)
 end if
 end subroutine complete_elliptic_integrals_standard_input_range
 
 elemental subroutine incomplete_elliptic_integrals_standard_input_range(phi,n,m,F,E,P,b,d,j)
 !!!!SJT: computes the incomplete elliptic integrals of kinds first to third given parameter m and characteristic n
 !!assumes 0<=phi<=pi/2, 0<=m<=1, 0<=n<=1
 use xelbdj2_all_routines, only: elbdj2
 
 !!input
 real(dp),intent(in) :: phi,n  !!elliptic amplitude, characteristic, and parameter
 real(qp),intent(in) :: m
 
 !!output
 real(dp),intent(out) :: F,E,P !!incomplete elliptic integrals of the first kind, second kind, and third kind
 real(dp),intent(out) :: b,d,j !!associated incomplete elliptic integrals
 
 !!internal variables
 real(dp) :: phic       !!complimentary Jacobi amplitude
 real(qp) :: mc         !!complimentary parameter
 real(qp) :: piio2_qp   !!pi/2 

 piio2_qp=pii_qp/2._qp 

 phic=real(piio2_qp-real(phi,qp),dp); ! use higher precision to ensure correctness of phi+phic=pi/2

 mc=1.0_qp-m
 call elbdj2(phi,phic,n,mc,b,d,j)  !!associated incomplete elliptic integrals
 F=b+d; E=b+real(mc,dp)*d; P=F+n*j !!build standard elliptic integrals from associated elliptic integrals
 end subroutine incomplete_elliptic_integrals_standard_input_range

 elemental subroutine incomplete_elliptic_integrals_standard_input_range_qp(phi,n,m,F,E,P,b,d,j)
 !!!!SJT: computes the incomplete elliptic integrals of kinds first to third given parameter m and characteristic n
 !!note: works in quadruple precision (accurate to about 25 digits)
 !!assumes 0<=phi<=pi/2, 0<=m<=1, 0<=n<=1
 use xelbdj2_all_routines, only: elbdj2_qp
 
 !!input
 real(qp),intent(in) :: phi,n  !!elliptic amplitude, characteristic, and parameter
 real(qp),intent(in) :: m
 
 !!output
 real(qp),intent(out) :: F,E,P !!incomplete elliptic integrals of the first kind, second kind, and third kind
 real(qp),intent(out) :: b,d,j !!associated incomplete elliptic integrals
 
 !!internal variables
 real(qp) :: piio2      !!constants
 real(qp) :: phic       !!complimentary variables
 real(qp) :: mc

 piio2=pii_qp/2._qp
 
 phic=piio2-phi; ! OG

 mc=1.0_qp-m
 call elbdj2_qp(phi,phic,n,mc,b,d,j)  !!associated incomplete elliptic integrals
 F=b+d; E=b+mc*d; P=F+n*j !!build standard elliptic integrals from associated elliptic integrals
 end subroutine incomplete_elliptic_integrals_standard_input_range_qp
 
 elemental subroutine incomplete_elliptic_integrals_standard_amp_large_parameter(phi,m,F,E)
 !!!!SJT: reduce parameter using reciprocal-modulus transformation (if needed)
 !!!!SJT: computes the incomplete elliptic integrals of first and second kind given amplitude phi and parameter m
 !!assumes 0<=phi<=pi/2 and 0<=m (m>1 is allowed)
 
 !!input
 real(dp),intent(in) :: phi,m !!elliptic amplitude and parameter
 
 !!output
 complex(dp),intent(out) :: F,E !!incomplete elliptic integrals of the first and second kind
 
 !!internal variables
 real(dp) :: piio2 !!constants
 !real(dp) :: mc    !!complimentary parameter
 real(dp) :: k !!elliptic modulus
 !real(dp) :: mrc !!compliment of the reciprocal parameter
 real(dp) :: Fc_temp,Ec_temp,Pc_temp !!temporary complete elliptic integrals
 real(dp) :: bc_temp,dc_temp,jc_temp !!temporary associated complete elliptic integrals
 real(dp) :: F_temp,E_temp,P_temp    !!temporary incomplete elliptic integrals
 real(dp) :: b_temp,d_temp,j_temp    !!temporary associated incomplete elliptic integrals
 real(qp) :: F_temp_qp,E_temp_qp,P_temp_qp    !!temporary incomplete elliptic integrals
 real(qp) :: b_temp_qp,d_temp_qp,j_temp_qp    !!temporary associated incomplete elliptic integrals
 real(dp) :: amp !!temporary elliptic amplitude variable
 real(dp) :: n !!characteristic
 real(qp) :: mr_qp,mrc_qp,m_qp,mc_qp !!variables related to elliptic parameter
 real(qp) :: u_qp!,amp_qp,sin_amp_qp !!argument of arcsine, elliptic amplitude, and sine amplitude 
 real(qp) :: k_qp

 real(qp) :: cos_theta,theta,sin_theta

 piio2=pii/2._dp
 n=0._dp !!arbitrary characteristic -- integrals of the third kind are not used
 
 if (phi.eq.piio2) then
  call complete_elliptic_integrals(m,F,E)
 else !!0<=phi<pi/2
  if (m.gt.1._dp) then !!large parameter m>1
   k_qp=sqrt(real(m,qp)); k=real(k_qp,dp)
   !mc=1.d0-m; mr=1.d0/m !!OG
   mc_qp=1.0_qp-real(m,qp)
   !mc=real(mc_qp,dp); mr_qp=1.0_qp/real(m,qp)
   mr_qp=1.0_qp/real(m,qp)
   !u=k*sin(phi) !!argument of arcsine
   !u_qp=real(k,qp)*sin(real(phi,qp))
   u_qp=k_qp*sin(real(phi,qp))
   if (u_qp.le.1.0_qp) then
    !amp=asin(u)
    amp=real(asin(u_qp),dp)
    call incomplete_elliptic_integrals_standard_input_range(amp,n,mr_qp,F_temp,E_temp,P_temp,b_temp,d_temp,j_temp)

    ! OG 
    !F=cmplx(F_temp/k,0.d0,dp); E=cmplx(k*E_temp+mc*F%RE,0.d0,dp)
    ! simplified Re[E] using associated incomplete elliptic integral definitions
    F=cmplx(F_temp/k,0._dp,dp); E=cmplx(b_temp/k,0._dp,dp)
   else
    !mrc=1.d0-mr; !!OG
    mrc_qp=1.0_qp-mr_qp; !mrc=real(mrc_qp,dp) !OG
    call complete_elliptic_integrals_standard_input_range(n,mr_qp,Fc_temp,Ec_temp,Pc_temp,bc_temp,dc_temp,jc_temp)

    ! OG    
    !sin_amp=(sqrt(u**2.d0-1.d0))/(u*sqrt(mrc)) ! OG
    !amp=asin(min(sin_amp,1.d0)) !!it is possible for sin_amp to be slightly greater than unity due to round off error
    !amp_qp=asin(min(sin_amp_qp,1.0_qp)); amp=real(amp_qp,dp)  ! best so far

    ! for trig forms (reduces numerical error)
    cos_theta=1._qp/tan(real(phi,qp))/sqrt(-mc_qp)
    theta=acos(cos_theta)
    sin_theta=sin(theta)

    ! OG
    !call incomplete_elliptic_integrals_standard_input_range(amp,n,mrc_qp,F_temp,E_temp,P_temp,b_temp,d_temp,j_temp)

    ! quadruple precision with trig forms to alleviate cancellation error below
    call incomplete_elliptic_integrals_standard_input_range_qp(theta,real(n,qp),mrc_qp,&
                                                              &F_temp_qp,E_temp_qp,P_temp_qp,&
                                                              &b_temp_qp,d_temp_qp,j_temp_qp)
    ! OG
    !F=cmplx(Fc_temp,-F_temp,dp)/k

    ! use quadruple precision integral in the imaginary part (eliminates siginicant error for large m and small |phi|)
    F=cmplx(Fc_temp,-F_temp_qp,dp)/k

    ! OG
    !E=cmplx(m*Ec_temp+mc*Fc_temp,&
    !  &-F_temp+m*E_temp+real(((1.0_qp-m)*sin_amp_qp*cos(amp_qp))/sqrt(1.0_qp-mrc_qp*sin_amp_qp**2),dp),dp)/k

    ! simplified Re[E] using associated incomplete elliptic integral definitions
    ! also added mc_qp to simplify notation
    !E=cmplx(bc_temp,&
    !  &-F_temp+m*E_temp+real((mc_qp*sin_amp_qp*cos(amp_qp))/sqrt(1.0_qp-mrc_qp*sin_amp_qp**2),dp),dp)/k

    ! simplified Im[E] using associated incomplete elliptic integral definitions
    !E=cmplx(bc_temp,&
    !  &-mc*b_temp+real((mc_qp*sin_amp_qp*cos(amp_qp))/sqrt(1.0_qp-mrc_qp*sin_amp_qp**2),dp),dp)/k

    ! applied trigonometric forms and added quad precision in Im[E] to reduce cancellation error
    E=cmplx(bc_temp,&
      &+mc_qp*(-b_temp_qp+sin_theta*cos_theta/sqrt(1._qp-mrc_qp*sin_theta**2_isp)),dp)/k
   end if
  else !!standard input ranges
   m_qp=real(m,qp)
   call incomplete_elliptic_integrals_standard_input_range(phi,n,m_qp,F_temp,E_temp,P_temp,b_temp,d_temp,j_temp)
   F=cmplx(F_temp,0._dp,dp); E=cmplx(E_temp,0._dp,dp)
  end if
 end if
 
 end subroutine incomplete_elliptic_integrals_standard_amp_large_parameter
 
 elemental subroutine incomplete_elliptic_integrals(phi,m,F,E)
 !!!!SJT: computes the incomplete elliptic integrals of first and second kind given amplitude phi and parameter m
 !!assumes phi is real
 !!assumes m>=0 (otherwise NaN is returned)
 use,intrinsic :: ieee_arithmetic,only: ieee_value,ieee_quiet_nan
 
 !!input
 real(dp),intent(in) :: phi,m !!elliptic amplitude and parameter
 
 !!output
 complex(dp),intent(out) :: F,E !!incomplete elliptic integrals of the first and second kind
 
 !!internal variables
 real(dp) :: piio2 !!constants
 complex(dp) :: Fc,Ec !!complete elliptic integral values computed for m>1
 real(dp) :: phi_std
 real(dp) :: N0
 !integer(isp) :: N,i_sign !!OG with single precision integer N
 integer(isp) :: i_sign
 integer(idp) :: N !! double precision integer N
 real(dp) :: nan !! NaN value

 ! validate m value (m<0 not currently supported)
 if (m.lt.0._dp) then ! return NaN values for unsupported m values
  nan=ieee_value(0._dp,ieee_quiet_nan)
  F=cmplx(nan,nan,dp); E=cmplx(nan,nan,dp)
  return
 end if

 piio2=pii/2._dp
 
 if ((0._dp.le.phi).and.(phi.le.piio2)) then !!standard amplitude range
  call incomplete_elliptic_integrals_standard_amp_large_parameter(phi,m,F,E)
 else !!amplitude outside of standard range
  N0=phi/pii
  !N=nint(N0,isp) !!integer multiplier -- possible overflow for large phi
  N=nint(N0,idp) !!integer multiplier -- double precision integer to handle large phi
  
  !phi_std=phi-N*pii !!find appropriate phi in standard range -- significant cancellation may occur when using double precision
  phi_std=real(real(phi,qp)-N*pii_qp,dp) !!find appropriate phi in standard range
  i_sign=1_isp
  if (phi_std.lt.0._dp) then
   i_sign=-1_isp
   phi_std=abs(phi_std)
  end if
 
  call incomplete_elliptic_integrals_standard_amp_large_parameter(phi_std,m,F,E)
  if (N.eq.0_idp) then ! slightly cheaper and correctly handles case for F where m --> 1
   F=i_sign*F
   E=i_sign*E
  else
   call complete_elliptic_integrals(m,Fc,Ec)
   F=2_isp*N*Fc+i_sign*F
   E=2_isp*N*Ec+i_sign*E
  end if
 end if
 
 end subroutine incomplete_elliptic_integrals

end module elliptic
