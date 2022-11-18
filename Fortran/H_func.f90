subroutine compute_H_func(x,z,t,lambda,k,zI,RaT,RaC,H_func)
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
real*8 :: H_func

!!internal variables
real*8 :: tol
real*8 :: xmid
real*8 :: xL,xR
complex*16 :: H,HL,HR

tol=1.d-8

xmid=lambda/2.d0

if ((x.le.0.d0).or.(x.ge.lambda).or.(z.le.0.d0).or.(z.ge.1.d0)) then !!if at or outside the domain boundaries return zero
 H_func=0.d0
elseif (((xmid-tol).le.x).and.(x.le.(xmid+tol))) then !!mid line
 xL=xmid-2.d0*tol; xR=xmid+2.d0*tol
 call compute_H(xL,z,t,lambda,k,zI,RaT,RaC,HL)
 call compute_H(xR,z,t,lambda,k,zI,RaT,RaC,HR)
 H_func=(HL%RE+HR%RE)/2.d0
else
 call compute_H(x,z,t,lambda,k,zI,RaT,RaC,H)
 H_func=H%RE
end if

end subroutine compute_H_func

subroutine compute_H(x,z,t,lambda,k,zI,RaT,RaC,H)
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
!write(*,*) "H t26=",t26
      t28 = (-x + lambda / 0.2D1)
      t29 = Heaviside(t28)
      t31 = 2 * t29 - 1
      t32 = (0, -1) * t31      
      t34 = t15 * t22                 !!D
      t35 = F(t)
      t36 = t35 * t11
      t37 = t36 * t34
      t39 = t37 * t1 * t32 + t26      !!eq. 27
!write(*,*) "H t37*t1*t32,t26=",t37*t1*t32,t26
      t40 = JacobiCN(t39, t25)
!write(*,*) "H cn t40=",t40
      t41 = 0.31415926535897932D1 / 2 + (0, 1) * log(sqrt(1 - t40 ** 2&
     &) + (0, 1) * t40)
      t46 = exp(-2 * (-t41 * t21 + zI) * k)
      t47 = 1 + t46
      t48 = t47 ** 2
      t49 = 0.1D1 / t48
      t54 = (0, 2) * t22 * t31 * k * 0.31415926535897932D1 * t49 * RaC
      t55 = t11 * t15
      t56 = dFdt(t)
      t58 = JacobiDN(t39, t25)
      t59 = JacobiSN(t39, t25)
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
!write(*,*) "H t92=",t92
      t93 = EllipticE(t25)
      t94 = t26 * t93
      t95 = EllipticK(t25)
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
!write(*,*) "H t148=",t148
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
