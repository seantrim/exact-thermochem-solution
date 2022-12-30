subroutine compute_H_func(x,z,t,lambda,k,zI,RaT,RaC,H_func)
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
real*8 :: H_func

!!internal variables
real*8 :: tol
real*8 :: xmid
real*8 :: xL,xR,xsym
complex*16 :: H,HL,HR

!!external functions
real*8 :: H_horizontal_boundaries 

tol=1.d-8

xmid=lambda/2.d0

if ((z.lt.0.d0).or.(z.gt.1.d0)) then !!beyond top or bottom boundaries
 H_func=0.d0
elseif ((z.eq.0.d0).or.(z.eq.1.d0)) then !!at top/bottom boundary
 H_func=H_horizontal_boundaries(z,RaC,RaT,k,zI)
elseif (x.lt.0.d0) then !!beyond left sidewall
 xsym=-x
 !call compute_H(xsym,z,t,lambda,k,zI,RaT,RaC,H)
 call compute_H_include_exterior(xsym,z,t,lambda,k,zI,RaT,RaC,H)
 H_func=H%RE
elseif (x.gt.lambda) then !!beyond right sidewall
 xsym=lambda-x
 !call compute_H(xsym,z,t,lambda,k,zI,RaT,RaC,H)
 call compute_H_include_exterior(xsym,z,t,lambda,k,zI,RaT,RaC,H)
 H_func=H%RE
elseif ((x.eq.0.d0).or.(x.eq.lambda)) then !!at sidewalls
 call compute_H_sidewalls(x,z,t,lambda,k,zI,RaT,RaC,H_func)
elseif (((xmid-tol).le.x).and.(x.le.(xmid+tol))) then !!mid line
 xL=xmid-2.d0*tol; xR=xmid+2.d0*tol
 !call compute_H(xL,z,t,lambda,k,zI,RaT,RaC,HL)
 !call compute_H(xR,z,t,lambda,k,zI,RaT,RaC,HR)
 call compute_H_include_exterior(xL,z,t,lambda,k,zI,RaT,RaC,HL)
 call compute_H_include_exterior(xR,z,t,lambda,k,zI,RaT,RaC,HR)
 H_func=(HL%RE+HR%RE)/2.d0
else
 !call compute_H(x,z,t,lambda,k,zI,RaT,RaC,H)
 call compute_H_include_exterior(x,z,t,lambda,k,zI,RaT,RaC,H)
 H_func=H%RE
end if

end subroutine compute_H_func

real*8 function H_horizontal_boundaries(z,RaC,RaT,k,zI)
implicit none

!!input
real*8 :: z
real*8 :: RaC,RaT
real*8 :: k,zI

!!internal variables
real*8 :: t1,t2,t3,t4,t5

t1=z-zI
t2=k*t1
t3=exp(2.d0*t2)
t4=t3*t3
t5=t3+1.d0

H_horizontal_boundaries=4.d0*(RaC*k**2*t3/t5**2 &
&- 2.d0*RaC*k**2*t4/t5**3)/RaT
end function H_horizontal_boundaries

subroutine compute_H_sidewalls(x,z,t,lambda,k,zI,RaT,RaC,H)
!!for H at the sidewalls
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
real*8 :: H

!!external functions
real*8 :: csc,cot,arccot
real*8 :: f_integral,f_func,f_derivative

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0
real*8 :: lam
real*8 :: H_ge,H_lt
real*8 :: Q,bigZ0
real*8 :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
real*8 :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33

lam=lambda

t9=pii*x/lam
t10=pii*z

t11=cos(t10)
t12=cos(t9)
t13=sin(t10)
t14=sin(t9)
t15=cot(t10)
t16=csc(t10)

t2=t15+t16
t1=log(abs(t2))

t17=f_func(t)
t18=f_integral(t)
t19=f_derivative(t)
t20=lam**2
t21=(t20 + 1)**2

t25=-pii*t12*t18/lam + t1
t3=exp(t25)
t4=1.d0/t3

t24=zI - arccot(-0.5d0*t4 + 0.5d0*t3)/pii
t22=k*t24
t23=k*(t24 - 1)

t5=exp(-2*t22)
t7=exp(-4*t22)
t6=exp(-2*t23)
t8=exp(-4*t23)

t26=(t15**2 + 1)
t27=(t4 - t3)**2
t28=(pii*t15*t16 + t26*pii)
t29=(pii**2*t15**2*t16 + 2*t26*pii**2*t15 + t26*pii**2*t16)
t30=(pii**2*t4*t18*t14/t20 + pii**2*t3*t18*t14/t20)**2*RaC*k
t31=(t5 + 1)
t32=(t6 + 1)
t33=(t27 + 4)

Q=t25

bigZ0=-(exp(-Q)-exp(Q))/2.d0
!write(*,*) "CHS z,bigZ0=",z,bigZ0

if (bigZ0.ge.0.d0) then
 H_ge=(t21*pii**4*t11*t12*t17/lam**3 - RaC + RaT + 4*RaC*k*(t28*t4/(t2) + t28*t3/(t2))*t5/(t33*pii*t31**2))*pii*t12*t17*t13/(RaT*&
     &lam) - (t21*pii**5*t12*t17*t13/lam**3 + 8*RaC*k*(t28*t4/(t2) + t28*t3/(t2))**2*(t4 - t3)*t5/(t33**2*pii*t31**2) - 4*RaC*k*(&
     &2*t28**2*t4/(t2)**2 - t29*t4/(t2) - t29*t3/(t2))*t5/(t33*pii*t31**2) - 16*RaC*k**2*(t28*t4/(t2) + t28*t3/(t2))**2*t5/(t33**&
     &2*pii**2*t31**2) + 32*RaC*k**2*(t28*t4/(t2) + t28*t3/(t2))**2*t7/(t33**2*pii**2*t31**3))/RaT - (t21*pii**5*t12*t17*t13/lam*&
     &*5 + 8*t30*(t4 - t3)*t5/(t33**2*pii*t31**2) - 16*t30**2*t5/(t33**2*pii**2*t31**2) - 4*(pii**4*t4*t18**2*t14**2/lam**4 - pii&
     &**4*t3*t18**2*t14**2/lam**4 - pii**3*t12*t4*t18/lam**3 - pii**3*t12*t3*t18/lam**3)*RaC*k*t5/(t33*pii*t31**2) + 32*t30**2*t7&
     &/(t33**2*pii**2*t31**3))/RaT - (t21*pii**3*t12*t13*t19/lam**3 + 4*(pii*t12*t4*t17/lam + pii*t12*t3*t17/lam)*RaC*k*t5/(t33*p&
     &ii*t31**2))/RaT
 H=H_ge
elseif (bigZ0.lt.0.d0) then
 H_lt=(t21*pii**4*t11*t12*t17/lam**3 - RaC + RaT + 4*RaC*k*(t28*t4/(t2) + t28*t3/(t2))*t6/(t33*pii*t32**2))*pii*t12*t17*t13/(RaT*&
     &lam) - (t21*pii**5*t12*t17*t13/lam**3 + 8*RaC*k*(t28*t4/(t2) + t28*t3/(t2))**2*(t4 - t3)*t6/(t33**2*pii*t32**2) - 4*RaC*k*(&
     &2*t28**2*t4/(t2)**2 - t29*t4/(t2) - t29*t3/(t2))*t6/(t33*pii*t32**2) - 16*RaC*k**2*(t28*t4/(t2) + t28*t3/(t2))**2*t6/(t33**&
     &2*pii**2*t32**2) + 32*RaC*k**2*(t28*t4/(t2) + t28*t3/(t2))**2*t8/(t33**2*pii**2*t32**3))/RaT - (t21*pii**5*t12*t17*t13/lam*&
     &*5 + 8*t30*(t4 - t3)*t6/(t33**2*pii*t32**2) - 16*t30**2*t6/(t33**2*pii**2*t32**2) - 4*(pii**4*t4*t18**2*t14**2/lam**4 - pii&
     &**4*t3*t18**2*t14**2/lam**4 - pii**3*t12*t4*t18/lam**3 - pii**3*t12*t3*t18/lam**3)*RaC*k*t6/(t33*pii*t32**2) + 32*t30**2*t8&
     &/(t33**2*pii**2*t32**3))/RaT - (t21*pii**3*t12*t13*t19/lam**3 + 4*(pii*t12*t4*t17/lam + pii*t12*t3*t17/lam)*RaC*k*t6/(t33*p&
     &ii*t32**2))/RaT
 H=H_lt
end if
end subroutine compute_H_sidewalls

subroutine compute_H_include_exterior(x,z,t,lambda,k,zI,RaT,RaC,H)
!!Based on output from Maple's CodeGeneration[Fortran] command using optimize=true (not optimize=tryhard)
!!Edits mad to Maple's output to allow routine to be compiled
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
complex*16 :: H

!!external functions
real*8 :: Heaviside,Dirac,Dirac_derivative,signum
real*8 :: F,dFdt,d2Fdt2,d3Fdt3
complex*16 :: JacobiSN,JacobiCN,JacobiDN
complex*16 :: InverseJacobiAM,EllipticK,EllipticE,JacobiZeta

!!internal functions
complex*16 :: t1,t2,t3,t5,t8,t11,t12,t13,t14,t15
complex*16 :: t16,t17,t21,t22,t23,t24,t25,t26,t28,t29,t31,t32,t34,t35,t36,t38,t39,t40,t45,t46,t47,t48,t53,t55,t56,t57,t58,t59,t60
complex*16 :: t62,t67,t69,t70,t72,t73,t74,t79,t80,t81,t83,t84,t85,t86,t87,t88,t89,t90,t91,t93,t95,t96,t97,t102,t103,t104,t105
complex*16 :: t106,t107,t108,t109,t112,t113,t114,t117,t121,t122,t124,t125,t127,t128,t130,t131,t132,t133,t134,t135,t136,t138,t140
complex*16 :: t141,t142,t144,t145,t147,t148,t151,t160,t167,t168,t169,t170,t173,t176,t179,t180,t182,t183,t186,t196,t200,t205,t206
complex*16 :: t207,t208,t209,t210,t211,t212,t216,t217,t220,t221,t222,t228,t232,t237,t239,t240,t241,t246,t247,t248,t254,t256,t264
complex*16 :: t266,t272,t273,t284,t285,t289,t291,t297,t300,t302,t303,t304,t310,t315,t317,t321,t335,t342,t355,t357,t365,t369,t374
complex*16 :: t378,t382,t386,t393,t399,t404,t407,t408,t410,t415,t417,t419,t420,t431,t436,t437,t452,t459,t461,t470,t480,t482,t484
complex*16 :: t488,t489,t497,t498,t505,t506,t510,t513,t516,t517,t526,t527,t531,t534,t536,t547,t554,t562,t590,t597,t617,t619,t652

      t1 = (0.31415926535897932D1 ** 2)
      t2 = 0.31415926535897932D1 * t1
      t3 = (lambda ** 2)
      t5 = (t3 + 1) ** 2
      t8 = 0.1D1 / lambda / t3
      t11 = 0.1D1 / lambda
      t12 = t11 * x * 0.31415926535897932D1
      t13 = cos(t12)
      t14 = (0.31415926535897932D1 * z)
      t15 = sin(t14)
      t16 = t15 * t13
      t17 = d2Fdt2(t) !diff(diff(F(t), t), t)
      t21 = 0.1D1 / 0.31415926535897932D1
      t22 = sin(t12)
      t23 = t15 * t22
      t24 = abs(t23)  !!D function
      t25 = 0.1D1 / t24
      t26 = InverseJacobiAM(t14, t25)
      t28 = (-x + lambda / 0.2D1)
      t29 = Heaviside(t28)
      t31 = 2 * t29 - 1
      t32 = (0, -1) * t31
      t34 = t11 * t24
      t35 = F(t)
      t36 = t35 * t34
      t38 = t36 * t1 * t32 + t26
      t39 = JacobiCN(t38, t25)
      t40 = 0.31415926535897932D1 / 2 + (0, 1) * log(sqrt(1 - t39 ** 2&
     &) + (0, 1) * t39)
      t45 = exp(-2 * (-t40 * t21 + zI) * k)
      t46 = 1 + t45
      t47 = t46 ** 2
      t48 = 0.1D1 / t47
      t53 = dFdt(t) !diff(F(t), t)
      t55 = JacobiDN(t38, t25)
      t56 = JacobiSN(t38, t25)
      t57 = t56 * t55
      t58 = t39 ** 2
      t59 = -t58 + 1
      t60 = sqrt(t59)
      t62 = t45 / t60
      t67 = 0.1D1 / RaT
      t69 = 0.31415926535897932D1 * t22
      t70 = cos(t14)
      t72 = t1 ** 2
      t73 = t5 * t72
      t74 = t3 ** 2
      t79 = t48 * RaC
      t80 = k * t79
      t81 = t24 ** 2
      t83 = 0.1D1 / t24 / t81
      t84 = 0.31415926535897932D1 * t83
      t85 = t11 * t84
      t86 = signum(t23) !abs(1, t23)
      t87 = JacobiZeta(t26, t25)
      t88 = EllipticE(t25)
      t89 = t26 * t88
      t90 = EllipticK(t25)
      t91 = 0.1D1 / t90
      t93 = t91 * t89 + t87
      t95 = 0.1D1 / t81
      t96 = -1 + t95
      t97 = 0.1D1 / t96
      t102 = 2 * t14
      t103 = sin(t102)
      t104 = t97 * t103
      t105 = t15 ** 2
      t106 = t105 * t95
      t107 = 1 - t106
      t108 = sqrt(t107)
      t109 = 0.1D1 / t108
      t112 = -t97 * t81 * t93 / 2 - t81 * t26 / 2 + t109 * t104 / 4
      t113 = t112 * t86
      t114 = t113 * t16
      t117 = Dirac(t28)
      t121 = 0.1D1 / t3
      t122 = t121 * t2
      t124 = t35 * t86
      t125 = t124 * t16
      t127 = -2 * t114 * t85 + (0, 2) * t36 * t1 * t117 + t125 * t122 *& 
     &t32
      t128 = t55 * t127
      t130 = 0.31415926535897932D1 * t95
      t131 = t11 * t130
      t132 = t56 ** 2
      t133 = t25 * t132
      t134 = -t96
      t135 = 0.1D1 / t134
      t136 = t39 * t135
      t138 = t38 * t24
      t140 = t135 * t24
      t141 = JacobiZeta(t38, t25)
      t142 = t88 * t38
      t144 = t91 * t142 + t141
      t145 = t144 * t140
      t147 = -t136 * t133 - t138 * t57 + t145 * t57
      t148 = t147 * t86
      t151 = -t148 * t16 * t131 - t56 * t128
      t160 = t11 * 0.31415926535897932D1
      t167 = t109 * 0.31415926535897932D1
      t168 = t22 * t83
      t169 = 0.31415926535897932D1 * t168
      t170 = t86 * t70
      t173 = 2 * t112 * t170 * t169
      t176 = t35 * t11
      t179 = t176 * t170 * t22 * t2 * t32 + t167 - t173
      t180 = t55 * t179
      t182 = t22 * t95
      t183 = 0.31415926535897932D1 * t182
      t186 = -t147 * t170 * t183 - t56 * t180
      t196 = t5 * 0.31415926535897932D1 * t72
      t200 = t53 * t16
      t205 = (k ** 2)
      t206 = t205 / t46 / t47 * RaC
      t207 = 0.1D1 / t1
      t208 = t151 ** 2
      t209 = t208 * t207
      t210 = 0.1D1 / t59
      t211 = t45 ** 2
      t212 = t211 * t210
      t216 = t81 ** 2
      t217 = 0.1D1 / t216
      t220 = t13 ** 2
      t221 = t105 * t220
      t222 = t86 ** 2
      t228 = t121 * t1 * t83
      t232 = 2*Dirac(t23) !signum(1, t23)
      t237 = t13 * t11
      t239 = t86 * t15
      t240 = t91 * t88
      t241 = 1 - t106 - t240
      t246 = 0.31415926535897932D1 * t25
      t247 = t237 * t246
      t248 = t88 ** 2
      t254 = t90 ** 2
      t256 = t95 * t15
      t264 = -t248 * t26 + t88 * (-t95 + 2 - t106) * t90 * t26 + (-t256& 
     &* t108 * t70 + t87 * t96 + (t26 * t96 + t87) * t107) * t254
      t266 = 0.1D1 / t254
      t272 = t24 * t90
      t273 = t24 * t88 - t272
      t284 = t95 * t266
      t285 = 0.31415926535897932D1 * t284
      t289 = t24 * t135 * t88 - t272
      t291 = t289 * t239 * t237
      t297 = t24 * t93
      t300 = t239 * t237
      t302 = t25 * t93
      t303 = t96 ** 2
      t304 = 0.1D1 / t303
      t310 = t24 * t26
      t315 = t83 * t109 * t304 * t103
      t317 = t86 * t16 * t160
      t321 = 0.1D1 / t108 / t107
      t335 = Dirac_derivative(t28) !Dirac(1, t28)
      t342 = t72 * t31
      t355 = t56 * t39
      t357 = t39 * t25
      t365 = t144 * t135 * t56 * t357 - t135 * t55 * t133 - t38 * t56 *& 
     &t357
      t369 = -t365 * t86 * t16 * t131 - t355 * t95 * t127
      t374 = t135 * t25
      t378 = t55 * t39
      t382 = -t56 * t132 * t374 - t144 * t378 * t140 + t378 * t138 + t56&
     & * t374
      t386 = -t382 * t86 * t16 * t131 + t55 * t39 * t127
      t393 = t121 * t1 * t95
      t399 = t25 * t56
      t404 = t136 * t95 * t132
      t407 = t134 ** 2
      t408 = 0.1D1 / t407
      t410 = t39 * t408 * t217 * t132
      t415 = t56 * t369
      t417 = t386 * t55
      t419 = t160 * t57
      t420 = t38 * t86
      t431 = t408 * t95
      t436 = t55 ** 2
      t437 = t436 - t240
      t452 = -t248 * t38 + t88 * (-t95 + t436 + 1) * t90 * t38 + (-t95 *&
     & t56 * t378 + t141 * t96 + (t38 * t96 + t141) * t436) * t254
      t459 = t95 * t38
      t461 = t273 * t86
      t470 = -2 * t386 * t136 * t399 + t317 * t404 + 2 * t317 * t410 - t&
     &151 * t135 * t133 - t138 * t415 - t138 * t417 - t420 * t16 * t419 &
     &- t127 * t24 * t57 + t145 * t415 + t145 * t417 + t144 * t135 * t86&
     & * t16 * t419 - 2 * t317 * t144 * t431 * t57 + (-t91 * t461 * t16 &
     &* t160 * t459 - t266 * t97 * t452 * t239 * t247 + t91 * t88 * t127&
     & + t291 * t285 * t142 + t437 * t127) * t140 * t57
      t480 = t21 * k * t79
      t482 = 0.1D1 / t60 / t59
      t484 = t39 * t45
      t488 = t205 * t79
      t489 = t45 * t210
      t497 = t186 ** 2
      t498 = t497 * t207
      t505 = t70 * 0.31415926535897932D1
      t506 = t86 * t505
      t510 = 2 * t506 * t22 * t105 * t83 - 2 * t505 * t256
      t513 = t22 ** 2
      t516 = t70 ** 2
      t517 = t222 * t516
      t526 = t1 * t513 * t83
      t527 = t232 * t516
      t531 = t167 - t173
      t534 = t505 * t22 * t25
      t536 = t266 * t97
      t547 = t289 * t170 * t69
      t554 = t170 * t69
      t562 = cos(t102)
      t590 = -t365 * t170 * t183 - t355 * t95 * t179
      t597 = -t382 * t170 * t183 + t55 * t39 * t179
      t617 = t56 * t590
      t619 = t597 * t55
      t652 = -2 * t597 * t136 * t399 + t554 * t404 + 2 * t554 * t410 - t&
     &186 * t135 * t133 - t138 * t617 - t138 * t619 - t420 * t505 * t22 &
     &* t57 - t179 * t24 * t57 + t145 * t617 + t145 * t619 + t144 * t135&
     & * t170 * t69 * t57 - 2 * t506 * t22 * t144 * t431 * t57 + (-t91 *&
     & t273 * t170 * t69 * t459 - t536 * t452 * t86 * t534 + t547 * t284&
     & * t142 + t91 * t88 * t179 + t437 * t179) * t140 * t57
      H = t67 * (-t17 * t16 * t8 * t5 * t2 + (0, 2) * t62 * t57 * t53 * &
     &t34 * t31 * 0.31415926535897932D1 * k * t48 * RaC) + t67 * (t53 * &
     &t23 / t74 * t73 + 2 * t62 * t151 * t21 * t80) * t53 * t70 * t69 - &
     &t67 * (-t53 * t70 * t13 * t8 * t73 + 2 * t62 * t186 * t21 * t80 + &
     &RaC - RaT) * t53 * t15 * t13 * t160 - t67 * (t200 / lambda / t74 *&
     & t196 + 8 * t212 * t209 * t206 + 2 * t62 * (-t56 * t55 * (6 * t112&
     & * t222 * t221 * t121 * t1 * t217 + 2 * t113 * t23 * t228 - 2 * t1&
     &12 * t232 * t221 * t228 - 2 * (-t97 * t81 * (-2 * t91 * t113 * t16&
     & * t160 * t83 * t88 - t91 * t26 * t273 * t239 * t237 * t130 - 2 * &
     &t241 * t112 * t239 * t237 * t84 - t266 * t97 * t264 * t239 * t247 &
     &+ t291 * t285 * t89) / 2 - t300 * 0.31415926535897932D1 * t97 * t2&
     &97 - t300 * 0.31415926535897932D1 * t304 * t302 + t114 * t11 * t24&
     &6 - t300 * 0.31415926535897932D1 * t310 + t317 * t315 / 2 - t86 * &
     &t237 * 0.31415926535897932D1 * t15 * t105 * t83 * t321 * t104 / 4)&
     & * t86 * t16 * t85 + (0, -2) * t36 * t1 * t335 + (0, 4) * t125 * t&
     &122 * t117 + (0, 1) * t124 * t23 * t8 * t342 + t35 * t232 * t221 *&
     & t8 * t72 * t32) - t56 * t369 * t127 - t386 * t128 + 2 * t147 * t2&
     &22 * t221 * t228 + t148 * t23 * t393 - t147 * t232 * t221 * t393 -&
     & t470 * t86 * t16 * t131) * t21 * t80 + 2 * t484 * t482 * t208 * t&
     &480 - 4 * t489 * t209 * t488) - t67 * (t200 * t8 * t196 + 8 * t212&
     & * t498 * t206 + 2 * t62 * (-t56 * t55 * (-t510 * t321 * 0.3141592&
     &6535897932D1 / 2 + 6 * t112 * t517 * t1 * t513 * t217 + 2 * t112 *&
     & t239 * t1 * t168 - 2 * t112 * t527 * t526 - 2 * (-t97 * t81 * (-t&
     &91 * t26 * t461 * t505 * t182 - t536 * t264 * t86 * t534 + t547 * &
     &t284 * t89 + t91 * t531 * t88 + t241 * t531) / 2 - t554 * t97 * t2&
     &97 - t554 * t304 * t302 - t81 * t531 / 2 - t506 * t22 * t310 + t10&
     &9 * t97 * t562 * 0.31415926535897932D1 / 2 + t554 * t315 / 2 - t51&
     &0 * t321 * t104 / 8) * t170 * t169 + (0, 1) * t176 * t239 * t22 * &
     &t342 + t176 * t527 * t513 * t72 * t32) - t56 * t590 * t179 - t597 &
     &* t180 + 2 * t147 * t517 * t526 + t147 * t239 * t1 * t182 - t147 *&
     & t527 * t1 * t513 * t95 - t652 * t170 * t183) * t21 * t80 + 2 * t4&
     &84 * t482 * t497 * t480 - 4 * t489 * t498 * t488)
end subroutine compute_H_include_exterior

subroutine compute_H(x,z,t,lambda,k,zI,RaT,RaC,H) !!OG
!!Based on output from Maple's CodeGeneration[Fortran] command using optimize=true (not optimize=tryhard)
!!Edits mad to Maple's output to allow routine to be compiled
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
complex*16 :: H

!!external functions
real*8 :: Heaviside,Dirac,Dirac_derivative
real*8 :: F,dFdt,d2Fdt2,d3Fdt3
complex*16 :: JacobiSN,JacobiCN,JacobiDN
complex*16 :: InverseJacobiAM,EllipticK,EllipticE,JacobiZeta

!!internal functions
complex*16 :: t1,t2,t3,t5,t8,t9,t11,t12,t13,t14,t15,t16,t17,t21,t22,t23,t24
complex*16 :: t25,t26,t28,t29,t31,t32,t34,t35,t36,t37,t39,t40,t41,t46,t47,t48,t49
complex*16 :: t54,t55,t56,t58,t59,t60,t61,t62,t63,t64,t65,t66,t70,t73,t74,t75
complex*16 :: t76,t77,t82,t83,t84,t86,t87,t88,t90,t91,t92,t93,t94,t95,t96,t98
complex*16 :: t99,t100,t101,t102,t103,t104,t107,t110,t111,t112,t113,t114,t115,t118
complex*16 :: t119,t122,t126,t127,t130,t132,t133,t135,t136,t137,t138,t139,t140,t141
complex*16 :: t144,t146,t147,t148,t149,t151,t152,t154,t157,t170,t172,t173,t174,t177
complex*16 :: t181,t182,t184,t187,t196,t202,t203,t204,t206,t209,t210,t211,t212,t213
complex*16 :: t214,t215,t226,t233,t246,t248,t253,t255,t263,t267,t269,t270,t271,t275
complex*16 :: t276,t279,t280,t284,t288,t289,t294,t295,t296,t302,t311,t313,t314,t317
complex*16 :: t322,t323,t325,t329,t335,t337,t340,t342,t350,t354,t355,t357,t363,t380
complex*16 :: t387,t396,t405,t408,t411,t425,t428,t431,t441,t442,t449,t452,t453,t460
complex*16 :: t462,t475,t479,t483,t499,t507,t515,t524,t536,t537,t541,t542,t544,t545
complex*16 :: t552,t554,t564,t567,t586,t609,t616,t638,t640,t671

      t1 = (0.31415926535897932D1 ** 2)
      t2 = 0.31415926535897932D1 * t1
      t3 = (lambda ** 2)
      t5 = (t3 + 1) ** 2
      t8 = 0.1D1 / lambda / t3
      t9 = t8 * t5 * t2
      t11 = 0.1D1 / lambda
      t12 = t11 * x * 0.31415926535897932D1
      t13 = cos(t12)
      t14 = (0.31415926535897932D1 * z)
      t15 = sin(t14)
      t16 = t15 * t13
      t17 = d2Fdt2(t)
      t21 = 0.1D1 / 0.31415926535897932D1
      t22 = sin(t12)
      t23 = 0.1D1 / t22
      t24 = 0.1D1 / t15
      t25 = t24 * t23                 !!1/D
!write(*,*) "H t14,t25=",t14,t25
      t26 = InverseJacobiAM(t14, t25) !!F(pi*z|1/D^2)
!write(*,*) "H IJAM t26=",t26
      t28 = (-x + lambda / 0.2D1)
!write(*,*) "H t28=",t28
      t29 = Heaviside(t28)
!write(*,*) "H Heaviside t29=",t29
      t31 = 2 * t29 - 1
      t32 = (0, -1) * t31      
      t34 = t15 * t22                 !!D
      t35 = F(t)
!write(*,*) "H F(t)=",F(t)
!write(*,*) "H dFdt(t)=",dFdt(t)
!write(*,*) "H d2Fdt2(t)=",d2Fdt2(t)
!write(*,*) "H d3Fdt3(t)=",d3Fdt3(t)
      t36 = t35 * t11
      t37 = t36 * t34
      t39 = t37 * t1 * t32 + t26      !!eq. 27
!write(*,*) "H t37*t1*t32,t26=",t37*t1*t32,t26
      call compute_JacobiSN_CN_DN(t39,t25,t59,t40,t58) !!SJT: avoid repeated computations
      !t40 = JacobiCN(t39, t25) !!OG
      t41 = 0.31415926535897932D1 / 2 + (0, 1) * log(sqrt(1 - t40 ** 2&
     &) + (0, 1) * t40)
      t46 = exp(-2 * (-t41 * t21 + zI) * k)
      t47 = 1 + t46
      t48 = t47 ** 2
      t49 = 0.1D1 / t48
      t54 = (0, 2) * t22 * t31 * k * 0.31415926535897932D1 * t49 * RaC
      t55 = t11 * t15
      t56 = dFdt(t)
      !t58 = JacobiDN(t39, t25) !!OG
      !t59 = JacobiSN(t39, t25) !!OG
!write(*,*) "H t39,t25=",t39,t25
!write(*,*) "H sn t59=",t59
!write(*,*) "H cn t40=",t40
!write(*,*) "H dn t58=",t58
      t60 = t59 * t58
      t61 = t40 ** 2
      t62 = -t61 + 1
      t63 = sqrt(t62)
      t64 = 0.1D1 / t63
      t65 = t46 * t64
      t66 = t65 * t60
      t70 = 0.1D1 / RaT
      t73 = cos(t14)
      t74 = t73 * 0.31415926535897932D1 * t22
      t75 = t1 ** 2
      t76 = t5 * t75
      t77 = t3 ** 2
      t82 = t49 * RaC
      t83 = k * t82
      t84 = t22 ** 2
      t86 = 0.1D1 / t22 / t84
      t87 = t15 ** 2
      t88 = 0.1D1 / t87
      t90 = 0.31415926535897932D1 * t88 * t86
      t91 = t13 * t11
!write(*,*) "H t26,t25=",t26,t25
      t92 = JacobiZeta(t26, t25)
!write(*,*) "H JZ t92=",t92
      call compute_EllipticK_EllipticE(t25,t95,t93) !!SJT: avoid repeated computations
      !t93 = EllipticE(t25) !!OG
      t94 = t26 * t93
      !t95 = EllipticK(t25) !!OG
!write(*,*) "H t25=",t25
!write(*,*) "H EK t95=",t95
!write(*,*) "H EE t93=",t93
      t96 = 0.1D1 / t95
      t98 = t96 * t94 + t92
      t99 = t84 * t98
      t100 = 0.1D1 / t84
      t101 = t88 * t100
      t102 = -1 + t101
      t103 = 0.1D1 / t102
      t104 = t103 * t87
      t107 = t84 * t26
      t110 = 2 * t14
      t111 = sin(t110)
      t112 = t103 * t111
      t113 = 1 - t100
      t114 = sqrt(t113)
      t115 = 0.1D1 / t114
      t118 = -t104 * t99 / 2 - t87 * t107 / 2 + t115 * t112 / 4
      t119 = t118 * t91
      t122 = Dirac(t28)
      t126 = t2 * t32
      t127 = 0.1D1 / t3
      t130 = t35 * t15 * t13 * t127
      t132 = -2 * t119 * t90 + (0, 2) * t37 * t1 * t122 + t130 * t126
      t133 = t58 * t132
      t135 = t24 * t100
      t136 = 0.31415926535897932D1 * t135
      t137 = t59 ** 2
      t138 = t23 * t137
      t139 = -t102
      t140 = 0.1D1 / t139
      t141 = t140 * t24
      t144 = t39 * t34
      t146 = t22 * t60
      t147 = t140 * t15
!write(*,*) "H t39,t25=",t39,t25
      t148 = JacobiZeta(t39, t25)
!write(*,*) "H  JZ t148=",t148
      t149 = t93 * t39
      t151 = t96 * t149 + t148
      t152 = t151 * t147
      t154 = -t40 * t141 * t138 - t144 * t60 + t152 * t146
      t157 = -t154 * t91 * t136 - t59 * t133
      t170 = t115 * 0.31415926535897932D1
      t172 = 0.1D1 / t15 / t87
      t173 = t172 * t100
      t174 = t73 * 0.31415926535897932D1
      t177 = 2 * t118 * t174 * t173
      t181 = t36 * t73 * t22 * t126 + t170 - t177
      t182 = t58 * t181
      t184 = t88 * t23
      t187 = -t154 * t174 * t184 - t59 * t182
      t196 = d3Fdt3(t)
      t202 = (k ** 2)
      t203 = t202 / t47 / t48 * RaC
      t204 = t31 ** 2
      t206 = t84 * t204 * t1
      t209 = t56 ** 2
      t210 = t209 * t127 * t87
      t211 = t58 ** 2
      t212 = t137 * t211
      t213 = 0.1D1 / t62
      t214 = t46 ** 2
      t215 = t214 * t213
      t226 = t209 * t127
      t233 = t204 * t2
      t246 = 0.1D1 / t63 / t62
      t248 = t40 * t46
      t253 = t202 * t82
      t255 = t46 * t213
      t263 = t5 * 0.31415926535897932D1 * t75
      t267 = t56 * t16
      t269 = 0.1D1 / t1
      t270 = t157 ** 2
      t271 = t270 * t269
      t275 = t84 ** 2
      t276 = 0.1D1 / t275
      t279 = t13 ** 2
      t280 = t279 * t127
      t284 = t127 * t1
      t288 = t96 * t93
      t289 = 1 - t100 - t288
      t294 = 0.31415926535897932D1 * t23
      t295 = t11 * t294
      t296 = t93 ** 2
      t302 = t95 ** 2
      t311 = -t296 * t26 + t93 * (-t101 + 2 - t100) * t95 * t26 + (-t135&
     & * t114 * t73 + t92 * t102 + (t26 * t102 + t92) * t113) * t302
      t313 = 0.1D1 / t302
      t314 = t313 * t103
      t317 = t11 * 0.31415926535897932D1
      t322 = t15 * t22 * t95
      t323 = t15 * t22 * t93 - t322
      t325 = t96 * t26
      t329 = 0.31415926535897932D1 * t88
      t335 = t100 * t313
      t337 = 0.31415926535897932D1 * t24
      t340 = t34 * t140 * t93 - t322
      t342 = t340 * t91 * t337
      t350 = 0.31415926535897932D1 * t103
      t354 = t102 ** 2
      t355 = 0.1D1 / t354
      t357 = t13 * t317
      t363 = t355 * t111
      t380 = Dirac_derivative(t28)
      t387 = t75 * t31
      t396 = t59 * t40 * t88
      t405 = t140 * t59
      t408 = -t140 * t58 * t137 * t25 + t151 * t405 * t40 * t25 - t39 *& 
     &t59 * t40 * t25
      t411 = -t396 * t100 * t132 - t408 * t91 * t136
      t425 = t58 * t40
      t428 = t58 * t40 * t15 * t22 * t39 - t59 * t137 * t140 * t25 - t15&
     &1 * t425 * t140 * t34 + t405 * t25
      t431 = t58 * t40 * t132 - t428 * t91 * t136
      t441 = t24 * t23 * t59
      t442 = t40 * t140
      t449 = t91 * 0.31415926535897932D1 * t40
      t452 = t139 ** 2
      t453 = 0.1D1 / t452
      t460 = t59 * t411
      t462 = t431 * t58
      t475 = t151 * t140
      t479 = t151 * t453
      t483 = t211 - t288
      t499 = -t296 * t39 + t93 * (-t101 + t211 + 1) * t95 * t39 + (-t88& 
     &* t100 * t59 * t425 + t148 * t102 + (t39 * t102 + t148) * t211) *& 
     &t302
      t507 = t96 * t323
      t515 = -2 * t431 * t442 * t441 + t449 * t141 * t100 * t137 + 2 * t&
     &449 * t453 * t172 * t276 * t137 - t157 * t141 * t138 - t144 * t460&
     & - t144 * t462 - t39 * t15 * t91 * 0.31415926535897932D1 * t60 - t&
     &132 * t34 * t60 + t152 * t22 * t460 + t152 * t22 * t462 + t475 * t&
     &16 * t317 * t60 - 2 * t357 * t479 * t135 * t60 + (-t507 * t91 * t3&
     &37 * t100 * t39 - t314 * t499 * t13 * t295 + t96 * t93 * t132 + t3&
     &42 * t335 * t149 + t483 * t132) * t147 * t146
      t524 = t21 * k * t82
      t536 = t187 ** 2
      t537 = t536 * t269
      t541 = t87 ** 2
      t542 = 0.1D1 / t541
      t544 = t73 ** 2
      t545 = t544 * t1
      t552 = t170 - t177
      t554 = t73 * t337
      t564 = t23 * t313
      t567 = t340 * t73 * t329
      t586 = cos(t110)
      t609 = -t396 * t100 * t181 - t408 * t174 * t184
      t616 = -t428 * t174 * t184 + t58 * t40 * t181
      t638 = t59 * t609
      t640 = t616 * t58
      t671 = -2 * t616 * t442 * t441 + t174 * t442 * t88 * t138 + 2 * t1&
     &74 * t40 * t453 * t542 * t86 * t137 - t187 * t141 * t138 - t144 * &
     &t638 - t144 * t640 - t39 * t174 * t146 - t181 * t34 * t60 + t152 *&
     & t22 * t638 + t152 * t22 * t640 + t475 * t174 * t146 - 2 * t174 * &
     &t479 * t184 * t60 + (-t507 * t174 * t88 * t23 * t39 - t313 * t103 &
     &* t499 * t554 + t567 * t564 * t149 + t96 * t93 * t181 + t483 * t18&
     &1) * t147 * t146
      H = t70 * (t66 * t56 * t55 * t54 - t17 * t16 * t9) + t70 * (t56 *& 
     &t34 / t77 * t76 + 2 * t65 * t157 * t21 * t83) * t56 * t74 - t70 * &
     &(-t56 * t73 * t13 * t8 * t76 + 2 * t65 * t187 * t21 * t83 + RaC - &
     &RaT) * t56 * t74 - t70 * (-t196 * t16 * t9 - 8 * t215 * t212 * t21&
     &0 * t206 * t203 + t66 * t17 * t55 * t54 - 2 * t46 * t64 * t137 * t&
     &40 * t226 * t204 * t2 * k * t82 + 2 * t65 * t40 * t211 * t210 * t8&
     &4 * t233 * t83 - 2 * t248 * t246 * t137 * t211 * t226 * t87 * t84 &
     &* t233 * t83 + 4 * t255 * t212 * t210 * t206 * t253) - t70 * (t267&
     & / lambda / t77 * t263 + 8 * t215 * t271 * t203 + 2 * t65 * (-t59 &
     &* t58 * (6 * t118 * t280 * t1 * t88 * t276 + 2 * t118 * t284 * t10&
     &1 - 2 * (-t104 * t84 * (-2 * t96 * t118 * t91 * t329 * t86 * t93 -&
     & t325 * t323 * t13 * t317 * t135 - 2 * t289 * t118 * t91 * t90 - t&
     &314 * t311 * t13 * t295 + t342 * t335 * t94) / 2 - t91 * t350 * t8&
     &7 * t22 * t98 - t357 * t355 * t23 * t98 + t119 * t294 - t357 * t87&
     & * t22 * t26 + t91 * t329 * t86 * t115 * t363 / 2 - t91 * 0.314159&
     &26535897932D1 * t86 / t114 / t113 * t112 / 4) * t91 * t90 + (0, -2&
     &) * t37 * t1 * t380 + (0, 4) * t130 * t2 * t122 + (0, 1) * t35 * t&
     &34 * t8 * t387) - t59 * t411 * t132 - t431 * t133 + 2 * t154 * t28&
     &0 * t1 * t24 * t86 + t154 * t284 * t25 - t515 * t91 * t136) * t21 &
     &* t83 + 2 * t248 * t246 * t270 * t524 - 4 * t255 * t271 * t253) - &
     &t70 * (t267 * t8 * t263 + 8 * t215 * t537 * t203 + 2 * t65 * (-t59&
     & * t58 * (6 * t118 * t545 * t542 * t100 + 2 * t118 * t1 * t101 - 2&
     & * (-t104 * t84 * (-t325 * t323 * t73 * 0.31415926535897932D1 * t1&
     &84 - t313 * t103 * t311 * t554 + t96 * t552 * t93 + t567 * t564 * &
     &t94 + t289 * t552) / 2 - t73 * t350 * t15 * t99 - t73 * 0.31415926&
     &535897932D1 * t355 * t24 * t98 - t87 * t84 * t552 / 2 - t73 * 0.31&
     &415926535897932D1 * t15 * t107 + t115 * t103 * t586 * 0.3141592653&
     &5897932D1 / 2 + t174 * t173 * t115 * t363 / 2) * t174 * t173 + (0,&
     & 1) * t35 * t55 * t22 * t387) - t59 * t609 * t181 - t616 * t182 + &
     &2 * t154 * t545 * t172 * t23 + t154 * t1 * t25 - t671 * t174 * t18&
     &4) * t21 * t83 + 2 * t248 * t246 * t536 * t524 - 4 * t255 * t537 *&
     & t253)

end subroutine compute_H
