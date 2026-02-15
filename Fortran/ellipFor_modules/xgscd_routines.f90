!!!!SJT: Update Notes
!!!! -This is an edited version of the file xgscd.txt 
!!!! -deleted original driver program to allow linkage of subroutines to other code
!!!! -deleted sample data at end of file so that the file could be compiled
!!!! -converted to free form -- comment character "!" used throughout -- "&" used at left and right sides for line continuation
!!!! -created the xgscd_routines module and moved procedures into the contains block to control outside access
!!!! -added use statement for kind_parameters module granting access to portable kind parameters in xgscd_routines module
!!!! -added implicit none statement to xgscd_routines module which extends to each routine in the contains block
!!!! -disabled save statement in variable declarations for thread safety
!!!! -added elemental keywords to all procedures (removed write statements that do not occur in practice)
!!!! -modified computations of m from mc (quad precision) in an attempt to reduce truncation error when mc is close to unity
!!!! -removed single precision routines
!!!! -variable declarations were modified to used the kind values from the kind_parameters module
!!!! -added intent(in) and intent(out) to argument declarations 
!!!! -replaced tab characters with spaces
!!!! -removed unused variables
!!!! -modernized syntax for parameter variable declarations
!!!! -in the code comments, "OG" is short for "Original code" and "SJT" indicates a modification
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xgscd_routines
 ! module for routines related to the computation of the pricipal Jacobi elliptic functions for standard input parameters
 use kind_parameters
 implicit none
 private
 public :: gscd ! computation of the principal Jacobi elliptic functions: sn, cn, and dn for standard input parameters
contains 
 elemental subroutine gscd(u,mc_qp,s,c,d) !!SJT
 !subroutine gscd(u,mc,s,c,d) !!SJT: original
 !
 !     Double precision subroutine to compute three Jacobian elliptic functions simultaneously
 !
 !     For general argument: -infty < u < infty
 !
 !     Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
 !       "Precise and Fast Computation of Jacobian Elliptic Functions by
 !        Conditional Duplication"
 !
 !     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
 !
 !     Used subprograms: scd2, elk
 !
 !     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
 !
 !     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
 !
 real(qp),intent(in) :: mc_qp !!SJT
 real(dp),intent(in) :: u !!SJT
 real(dp),intent(out) :: s,c,d !!SJT
 !real(dp) mc,m,kc,ux,k,kh,kh3,kh5,kh7,k2,k3,k4,sx !!SJT: OG
 real(dp) kc,ux,k,kh,kh3,kh5,kh7,k2,k3,k4,sx !!SJT: removed m and mc (mc_qp used instead)
 !real(dp) elk !!SJT: declaration of elk function not required due to xgscd_routines module 
 !
 !m=1.d0-mc !!SJT: original
 !m=real(1.0_qp-mc_qp,dp); !!SJT: removed 
 !mc=real(mc_qp,dp) !!SJT: quad precision used to reduce the impact of cancellation errors
 !kc=sqrt(mc) !! SJT: OG
 kc=real(sqrt(mc_qp),dp)
 ux=abs(u)
 if(ux.lt.0.785e0_dp) then
     !call scd2(ux,mc,s,c,d) !!SJT: OG
     call scd2(ux,mc_qp,s,c,d) !!SJT
 else
     !k=elk(mc) !!SJT: OG
     k=elk(mc_qp) !!SJT
     kh=k*0.5e0_dp; kh3=k*1.5e0_dp; kh5=k*2.5e0_dp; kh7=k*3.5e0_dp;
     k2=k*2.e0_dp; k3=k*3.e0_dp; k4=k*4.e0_dp
     !ux=ux-k4*dble(int(ux/k4)) !!SJT: OG
     ux=ux-k4*real(int(ux/k4,idp),dp) !!SJT: updated to prevent integer overflow
     if(ux.lt.kh) then
         !call scd2(ux,mc,s,c,d) !!SJT: OG
         call scd2(ux,mc_qp,s,c,d) !!SJT
     elseif(ux.lt.k) then
         ux=k-ux
         !call scd2(ux,mc,s,c,d) !!SJT: OG
         call scd2(ux,mc_qp,s,c,d) !!SJT
         sx=c/d; c=kc*s/d; s=sx; d=kc/d
     elseif(ux.lt.kh3) then
         ux=ux-k
         !call scd2(ux,mc,s,c,d) !!SJT: OG
         call scd2(ux,mc_qp,s,c,d) !!SJT
         sx=c/d; c=-kc*s/d; s=sx; d=kc/d
     elseif(ux.lt.k2) then
         ux=k2-ux
         !call scd2(ux,mc,s,c,d) !!SJT: OG
         call scd2(ux,mc_qp,s,c,d) !!SJT
         c=-c
     elseif(ux.lt.kh5) then
         ux=ux-k2
         !call scd2(ux,mc,s,c,d) !!SJT: OG
         call scd2(ux,mc_qp,s,c,d) !!SJT
         s=-s; c=-c
     elseif(ux.lt.k3) then
         ux=k3-ux
         !call scd2(ux,mc,s,c,d) !!SJT: OG
         call scd2(ux,mc_qp,s,c,d) !!SJT
         sx=-c/d; c=-kc*s/d; s=sx; d=kc/d
     elseif(ux.lt.kh7) then
         ux=ux-k3
         !call scd2(ux,mc,s,c,d) !!SJT: OG
         call scd2(ux,mc_qp,s,c,d) !!SJT
         sx=-c/d; c=kc*s/d; s=sx; d=kc/d
     else
         ux=k4-ux
         !call scd2(ux,mc,s,c,d) !!SJT: OG
         call scd2(ux,mc_qp,s,c,d) !!SJT
         s=-s
     endif
 endif
 if(u.lt.0.e0_dp) s=-s
 return
 end
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     elemental subroutine scd2(u,mc_qp,s,c,d) !!SJT
 !    subroutine scd2(u,mc,s,c,d) !!SJT: original
 !
 !        Double precision subroutine to compute three Jacobian elliptic functions simultaneously
 !
 !   For limited argument: 0 <= u < K/2
 !
 !     Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
 !       "Precise and Fast Computation of Jacobian Elliptic Functions by
 !        Conditional Duplication"
 !
 !     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
 !
 !     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
 !
 !     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
 !
     real(dp) :: x !!SJT
     real(dp),intent(in) :: u !!SJT
     real(dp),intent(out) :: s,c,d !!SJT
     real(dp) m,mc,uA,uT,u0,v,a,b,y,z,my,mc2,m2,xz,w !!SJT
     integer(isp) n,j,i
     real(dp),parameter :: B10=1.e0_dp/24.e0_dp,B11=1.e0_dp/6.e0_dp,B20=1.e0_dp/720.e0_dp !!SJT
     real(dp),parameter :: B21=11.e0_dp/180.e0_dp,B22=1.e0_dp/45.e0_dp !!SJT
     real(qp),intent(in) :: mc_qp !!SJT
     !m=1.d0-mc; uA=1.76269d0+mc*1.16357d0; uT=5.217d-3-m*2.143d-3; u0=u !!SJT: original
     m=real(1.0_qp-mc_qp,dp); mc=real(mc_qp,dp) !!SJT
     uA=1.76269e0_dp+mc*1.16357e0_dp; uT=5.217e-3_dp-m*2.143e-3_dp; u0=u !!SJT
     do n=0_isp,20_isp
         if(u0.lt.uT) goto 1
         u0=u0*0.5e0_dp
     enddo
     !write(*,*) "(scd2) Too large input argument: u=", u !!SJT does not occur in practice
 1 continue
     v=u0*u0; a=1.e0_dp; b=v*(0.5e0_dp-v*(B10+m*B11-v*(B20+m*(B21+m*B22))))
     if(u.lt.uA) then
         do j=1_isp,n
             y=b*(a*2.e0_dp-b); z=a*a; my=m*y; b=(y*2.e0_dp)*(z-my); a=z*z-my*y
         enddo
     else
         do j=1_isp,n
             y=b*(a*2.e0_dp-b); z=a*a; my=m*y
             if(z.lt.my*2.e0_dp) goto 2
             b=(y*2.e0_dp)*(z-my); a=z*z-my*y
         enddo
     endif
     b=b/a; y=b*(2.e0_dp-b); c=1.e0_dp-b; s=sqrt(y); d=sqrt(1.e0_dp-m*y)
      return
 2 continue
     c=a-b; mc2=mc*2.e0_dp; m2=m*2.e0_dp
     do i=j,n
         x=c*c; z=a*a; w=m*x*x-mc*z*z; xz=x*z; c=mc2*xz+w; a=m2*xz-w
     enddo
     c=c/a; x=c*c; s=sqrt(1.e0_dp-x); d=sqrt(mc+m*x)
     return; end
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       elemental real(dp) function elk(mc_qp) !!SJT
 !      real*8 function elk(mc) !!SJT: original
 !c
 !c        Double precision complete elliptic integral of the first kind
 !c
 !c     Reference: T. Fukushima, (2009) Celest. Mech. Dyn. Astron. 105, 305-328
 !c        "Fast Computation of Complete Elliptic Integrlals and Jacobian
 !c         Elliptic Functions"
 !c
 !c     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
 !c
 !c     Inputs: mc   = complementary parameter 0 <= mc   <= 1
 !c
 !c     Output: elk
 !c
         real(dp) mc
         real(dp) mcold,PIHALF,PIINV,elkold,TINY,m,mx
         real(dp) kkc,nome
         !real(dp) :: kc !!SJT
         real(qp),intent(in) :: mc_qp !!SJT
 !c
         real(dp),parameter :: D1=1.e0_dp/16.e0_dp,D2=1.e0_dp/32.e0_dp,D3=21.e0_dp/1024.e0_dp !!SJT  
         real(dp),parameter :: D4=31.e0_dp/2048.e0_dp,D5=6257.e0_dp/524288.e0_dp              !!SJT 
         real(dp),parameter :: D6=10293.e0_dp/1048576.e0_dp,D7=279025.e0_dp/33554432.e0_dp    !!SJT  
         real(dp),parameter :: D8=483127.e0_dp/67108864.e0_dp                                 !!SJT  
         real(dp),parameter :: D9=435506703.e0_dp/68719476736.e0_dp                           !!SJT 
         real(dp),parameter :: D10=776957575.e0_dp/137438953472.e0_dp                         !!SJT 
         real(dp),parameter :: D11=22417045555.e0_dp/4398046511104.e0_dp                      !!SJT 
         real(dp),parameter :: D12=40784671953.e0_dp/8796093022208.e0_dp                      !!SJT 
         real(dp),parameter :: D13=9569130097211.e0_dp/2251799813685248.e0_dp                 !!SJT 
         real(dp),parameter :: D14=17652604545791.e0_dp/4503599627370496.e0_dp                !!SJT 
 !c
         !logical first/.TRUE./ !!SJT: note that logical variable first is .true. at the start of each routine call
         !save first,mcold,PIHALF,PIINV,elkold,TINY !!SJT: disabled for thread safety
         logical first !!SJT: initialization cannot be in declaration statement for thread safety purposes
         first=.true. !!SJT: initialization outside of declaration for thread safety
 !c
         if(first) then
                 first=.FALSE.
                 mcold=1.e0_dp
                 PIHALF=atan(1.e0_dp)*2.e0_dp
                 PIINV=0.5e0_dp/PIHALF
                 elkold=PIHALF
                 TINY=1.e-99_dp
         endif
         !m=1.d0-mc !!SJT: original
         m=real(1.0_qp-mc_qp,dp); mc=real(mc_qp,dp) !!SJT
         !kc=sqrt(mc); m=(1.d0+kc)*(1.d0-kc) !!!!SJT: reduce truncation error
         if(abs(m).lt.1.e-16_dp) then
                 elk=PIHALF
         elseif(abs(mc-mcold).lt.1.11e-16_dp*mc) then
                 elk=elkold
         elseif(mc.lt.TINY) then
           elk=1.3862943611198906e0_dp-0.5e0_dp*log(TINY)
         elseif(mc.lt.1.11e-16_dp) then
           elk=1.3862943611198906e0_dp-0.5e0_dp*log(mc)
         elseif(mc.lt.0.1e0_dp) then
                 nome=mc*(D1+mc*(D2+mc*(D3+mc*(D4+mc*(D5+mc*(D6&
      &                +mc*(D7+mc*(D8+mc*(D9+mc*(D10+mc*(D11+mc*(D12&
      &                +mc*(D13+mc*D14)))))))))))))
                 mx=mc-0.05e0_dp
 !c
 !c        K'
 !c
                 kkc=1.591003453790792180e0_dp+mx*(&
      &                0.416000743991786912e0_dp+mx*(&
      &                0.245791514264103415e0_dp+mx*(&
      &                0.179481482914906162e0_dp+mx*(&
      &                0.144556057087555150e0_dp+mx*(&
      &                0.123200993312427711e0_dp+mx*(&
      &                0.108938811574293531e0_dp+mx*(&
      &                0.098853409871592910e0_dp+mx*(&
      &                0.091439629201749751e0_dp+mx*(&
      &                0.085842591595413900e0_dp+mx*(&
      &                0.081541118718303215e0_dp))))))))))
 !c
                 elk=-kkc*PIINV*log(nome)
         elseif(m.le.0.1e0_dp) then
                 mx=m-0.05e0_dp
                 elk=1.591003453790792180e0_dp+mx*(&
      &                0.416000743991786912e0_dp+mx*(&
      &                0.245791514264103415e0_dp+mx*(&
      &                0.179481482914906162e0_dp+mx*(&
      &                0.144556057087555150e0_dp+mx*(&
      &                0.123200993312427711e0_dp+mx*(&
      &                0.108938811574293531e0_dp+mx*(&
      &                0.098853409871592910e0_dp+mx*(&
      &                0.091439629201749751e0_dp+mx*(&
      &                0.085842591595413900e0_dp+mx*(&
      &                0.081541118718303215e0_dp))))))))))
         elseif(m.le.0.2e0_dp) then
                 mx=m-0.15e0_dp
                 elk=1.635256732264579992e0_dp+mx*(&
      &                0.471190626148732291e0_dp+mx*(&
      &                0.309728410831499587e0_dp+mx*(&
      &                0.252208311773135699e0_dp+mx*(&
      &                0.226725623219684650e0_dp+mx*(&
      &                0.215774446729585976e0_dp+mx*(&
      &                0.213108771877348910e0_dp+mx*(&
      &                0.216029124605188282e0_dp+mx*(&
      &                0.223255831633057896e0_dp+mx*(&
      &                0.234180501294209925e0_dp+mx*(&
      &                0.248557682972264071e0_dp+mx*(&
      &                0.266363809892617521e0_dp+mx*(&
      &                0.287728452156114668e0_dp))))))))))))
         elseif(m.le.0.3e0_dp) then
                 mx=m-0.25e0_dp
                 elk=1.685750354812596043e0_dp+mx*(&
      &                0.541731848613280329e0_dp+mx*(&
      &                0.401524438390690257e0_dp+mx*(&
      &                0.369642473420889090e0_dp+mx*(&
      &                0.376060715354583645e0_dp+mx*(&
      &                0.405235887085125919e0_dp+mx*(&
      &                0.453294381753999079e0_dp+mx*(&
      &                0.520518947651184205e0_dp+mx*(&
      &                0.609426039204995055e0_dp+mx*(&
      &                0.724263522282908870e0_dp+mx*(&
      &                0.871013847709812357e0_dp+mx*(&
      &                1.057652872753547036e0_dp)))))))))))
         elseif(m.le.0.4e0_dp) then
                 mx=m-0.35e0_dp
                 elk=1.744350597225613243e0_dp+mx*(&
      &                0.634864275371935304e0_dp+mx*(&
      &                0.539842564164445538e0_dp+mx*(&
      &                0.571892705193787391e0_dp+mx*(&
      &                0.670295136265406100e0_dp+mx*(&
      &                0.832586590010977199e0_dp+mx*(&
      &                1.073857448247933265e0_dp+mx*(&
      &                1.422091460675497751e0_dp+mx*(&
      &                1.920387183402304829e0_dp+mx*(&
      &                2.632552548331654201e0_dp+mx*(&
      &                3.652109747319039160e0_dp+mx*(&
      &                5.115867135558865806e0_dp+mx*(&
      &                7.224080007363877411e0_dp))))))))))))
         elseif(m.le.0.5e0_dp) then
                 mx=m-0.45e0_dp
                 elk=1.813883936816982644e0_dp+mx*(&
      &                0.763163245700557246e0_dp+mx*(&
      &                0.761928605321595831e0_dp+mx*(&
      &                0.951074653668427927e0_dp+mx*(&
      &                1.315180671703161215e0_dp+mx*(&
      &                1.928560693477410941e0_dp+mx*(&
      &                2.937509342531378755e0_dp+mx*(&
      &                4.594894405442878062e0_dp+mx*(&
      &                7.330071221881720772e0_dp+mx*(&
      &                11.87151259742530180e0_dp+mx*(&
      &                19.45851374822937738e0_dp+mx*(&
      &                32.20638657246426863e0_dp+mx*(&
      &                53.73749198700554656e0_dp+mx*(&
      &                90.27388602940998849e0_dp)))))))))))))
         elseif(m.le.0.6e0_dp) then
                 mx=m-0.55e0_dp
                 elk=1.898924910271553526e0_dp+mx*(&
      &                0.950521794618244435e0_dp+mx*(&
      &                1.151077589959015808e0_dp+mx*(&
      &                1.750239106986300540e0_dp+mx*(&
      &                2.952676812636875180e0_dp+mx*(&
      &                5.285800396121450889e0_dp+mx*(&
      &                9.832485716659979747e0_dp+mx*(&
      &                18.78714868327559562e0_dp+mx*(&
      &                36.61468615273698145e0_dp+mx*(&
      &                72.45292395127771801e0_dp+mx*(&
      &                145.1079577347069102e0_dp+mx*(&
      &                293.4786396308497026e0_dp+mx*(&
      &                598.3851815055010179e0_dp+mx*(&
      &                1228.420013075863451e0_dp+mx*(&
      &                2536.529755382764488e0_dp))))))))))))))
         elseif(m.le.0.7e0_dp) then
                 mx=m-0.65e0_dp
                 elk=2.007598398424376302e0_dp+mx*(&
      &                1.248457231212347337e0_dp+mx*(&
      &                1.926234657076479729e0_dp+mx*(&
      &                3.751289640087587680e0_dp+mx*(&
      &                8.119944554932045802e0_dp+mx*(&
      &                18.66572130873555361e0_dp+mx*(&
      &                44.60392484291437063e0_dp+mx*(&
      &                109.5092054309498377e0_dp+mx*(&
      &                274.2779548232413480e0_dp+mx*(&
      &                697.5598008606326163e0_dp+mx*(&
      &                1795.716014500247129e0_dp+mx*(&
      &                4668.381716790389910e0_dp+mx*(&
      &                12235.76246813664335e0_dp+mx*(&
      &                32290.17809718320818e0_dp+mx*(&
      &                85713.07608195964685e0_dp+mx*(&
      &                228672.1890493117096e0_dp+mx*(&
      &                612757.2711915852774e0_dp))))))))))))))))
         elseif(m.le.0.8e0_dp) then
                 mx=m-0.75e0_dp
                 elk=2.156515647499643235e0_dp+mx*(&
      &                1.791805641849463243e0_dp+mx*(&
      &                3.826751287465713147e0_dp+mx*(&
      &                10.38672468363797208e0_dp+mx*(&
      &                31.40331405468070290e0_dp+mx*(&
      &                100.9237039498695416e0_dp+mx*(&
      &                337.3268282632272897e0_dp+mx*(&
      &                1158.707930567827917e0_dp+mx*(&
      &                4060.990742193632092e0_dp+mx*(&
      &                14454.00184034344795e0_dp+mx*(&
      &                52076.66107599404803e0_dp+mx*(&
      &                189493.6591462156887e0_dp+mx*(&
      &                695184.5762413896145e0_dp+mx*(&
      &                2.567994048255284686e6_dp+mx*(&
      &                9.541921966748386322e6_dp+mx*(&
      &                3.563492744218076174e7_dp+mx*(&
      &                1.336692984612040871e8_dp+mx*(&
      &                5.033521866866284541e8_dp+mx*(&
      &                1.901975729538660119e9_dp+mx*(&
      &                7.208915015330103756e9_dp)))))))))))))))))))
         elseif(m.le.0.85e0_dp) then
                 mx=m-0.825e0_dp
                 elk=2.318122621712510589e0_dp+mx*(&
      &                2.616920150291232841e0_dp+mx*(&
      &                7.897935075731355823e0_dp+mx*(&
      &                30.50239715446672327e0_dp+mx*(&
      &                131.4869365523528456e0_dp+mx*(&
      &                602.9847637356491617e0_dp+mx*(&
      &                2877.024617809972641e0_dp+mx*(&
      &                14110.51991915180325e0_dp+mx*(&
      &                70621.44088156540229e0_dp+mx*(&
      &                358977.2665825309926e0_dp+mx*(&
      &                1.847238263723971684e6_dp+mx*(&
      &                9.600515416049214109e6_dp+mx*(&
      &                5.030767708502366879e7_dp+mx*(&
      &                2.654441886527127967e8_dp+mx*(&
      &                1.408862325028702687e9_dp+mx*(&
      &                7.515687935373774627e9_dp)))))))))))))))
       else
                 mx=m-0.875e0_dp
                 elk=2.473596173751343912e0_dp+mx*(&
      &                3.727624244118099310e0_dp+mx*(&
      &                15.60739303554930496e0_dp+mx*(&
      &                84.12850842805887747e0_dp+mx*(&
      &                506.9818197040613935e0_dp+mx*(&
      &                3252.277058145123644e0_dp+mx*(&
      &                21713.24241957434256e0_dp+mx*(&
      &                149037.0451890932766e0_dp+mx*(&
      &                1.043999331089990839e6_dp+mx*(&
      &                7.427974817042038995e6_dp+mx*(&
      &                5.350383967558661151e7_dp+mx*(&
      &                3.892498869948708474e8_dp+mx*(&
      &                2.855288351100810619e9_dp+mx*(&
      &                2.109007703876684053e10_dp+mx*(&
      &                1.566998339477902014e11_dp+mx*(&
      &                1.170222242422439893e12_dp+mx*(&
      &                8.777948323668937971e12_dp+mx*(&
      &                6.610124275248495041e13_dp+mx*(&
      &                4.994880537133887989e14_dp+mx*(&
      &                3.785974339724029920e15_dp)))))))))))))))))))
         endif
 !c
         mcold=mc
         elkold=elk
       return
       end
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module xgscd_routines
