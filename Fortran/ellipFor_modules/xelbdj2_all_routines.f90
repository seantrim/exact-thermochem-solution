!!!!SJT: Update Notes
!!!! -edited version of xelbdj2_all.txt for the purposes of compiling routines with other Fortran code
!!!! -removed original driver program
!!!! -removed sample data at end of file
!!!! -created the xelbdj2_all_routines module to facilitate access to the routines (now in the module contains block)
!!!! -explicit declarations for serj and uatan functions were removed because module definition already grants access
!!!! -added use statememt for kind_parameters module (in kind_parameters.f90) for access to portable kind parameters
!!!! -added implicit none statement in module definition, which extends to all procedures in the contains block
!!!! -modified computations of m from mc to reduce truncation error when mc is close to unity
!!!! -modified upper limits on two do loops in celbdj to avoid an index potentially going out of bounds
!!!! -disabled save and data statements in variable declarations for thread safety
!!!! -added elemental keywords to all procedures (removed write statements that do not occur in practice)
!!!! -variable declarations have been modified to use the kind parameters from the kind_parameters module
!!!! -added intent(in) and intent(out) to argument declarations
!!!! -replaced tab characters with spaces
!!!! -removed unused variables
!!!! -added quadruple precision versions of the original routines
!!!! -modernized syntax for parameter variable declarations
!!!! -in the code comments, "OG" is short for "Original code" and "SJT" indicates a modification
module xelbdj2_all_routines
 ! module containing subroutines for the calculation of associated Legendre elliptic integrals for standard input ranges
 use kind_parameters
 implicit none
 private
 public :: elbdj2 ! compute the associated Legendre integrals of the first, second, and third kinds for standard input ranges
 public :: elbdj2_qp ! extra precision version of elbdj2
contains
 !---------------------------------------------------------------------------
 elemental subroutine elbdj2(phi,phic,n,mc_qp,b,d,j) !!SJT
 !subroutine elbdj2(phi,phic,n,mc,b,d,j) !!SJT: OG
 !
 !     Simultaneous computation of associate elliptic integrals
 !     of third kind, B(phi|m), D(phi|m), and J(phi,n|m)
 !
 !     Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !       Precise and fast computation of a general incomplete elliptic integral
 !       of third kind by half and double argument transformations
 !
 !     Modified Version by employing "elsbdj" only
 !
 !     Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
 !
 !     Used subprograms: celbd,celbdj,elsbdj,serbd,serj,uatan
 !
 !     Inputs: phi  = argument                0 <= phi  <= PI/2
 !             phic = complementar argument   0 <= phic <= PI/2
 !             n    = characteristic          0 <= n    <= 1
 !             mc   = complementary parameter 0 <= mc   <= 1
 !
 !     Outputs: b, d, j
 !
 !     CAUTION: phi and phic must satisfy condition, phi + phic = PI/2
 !
 real(dp),intent(in) :: phi,phic,n !!SJT
 real(dp),intent(out) :: b,d,j !!SJT
 real(dp) mc,m,nc,h,c,x,d2,z,bc,dc,jc,sz,t,v,t2 !!SJT
 !real(dp) uatan !!SJT: due to module definition, this routine already has access to uatan function 
 real(qp),intent(in) :: mc_qp !!SJT
 !
 if(phi.lt.1.249e0_dp) then   ! Modified: Pascal Leroy's suggestion 2019/08/12
     !call elsbdj(sin(phi),n,mc,b,d,j) !!SJT: OG
     call elsbdj(sin(phi),n,mc_qp,b,d,j) !!SJT
 else
     !m=1.d0-mc !!SJT: original
     m=real(1.0_qp-mc_qp,dp); mc=real(mc_qp,dp) !!SJT
     nc=1.e0_dp-n
     h=n*nc*(n-m)
     c=sin(phic)
     x=c*c
     d2=mc+m*x
     if(x.lt.0.9e0_dp*d2) then
         z=c/sqrt(d2)
         !call elsbdj(z,n,mc,b,d,j) !!SJT: OG
         call elsbdj(z,n,mc_qp,b,d,j) !!SJT
         !call celbdj(nc,mc,bc,dc,jc) !!SJT: original
         call celbdj(nc,mc_qp,bc,dc,jc) !!SJT
                 sz=z*sqrt(1.e0_dp-x)
                 t=sz/nc
                 b=bc-(b-sz)
                 d=dc-(d+sz)
                 j=jc-(j+uatan(t,h))
         else
                 v=mc*(1.e0_dp-x)
                 if(v.lt.x*d2) then
             !call elsbdj(sqrt(1.e0_dp-c*c),n,mc,b,d,j) !!SJT: OG
             call elsbdj(sqrt(1.e0_dp-c*c),n,mc_qp,b,d,j) !!SJT
         else
             t2=(1.e0_dp-x)/d2
             !call elsbdj(sqrt(1.e0_dp-mc*t2),n,mc,b,d,j) !!SJT: OG
             call elsbdj(sqrt(1.e0_dp-mc*t2),n,mc_qp,b,d,j) !!SJT
             !call celbdj(nc,mc,bc,dc,jc) !!SJT: OG
             call celbdj(nc,mc_qp,bc,dc,jc) !!SJT
                     sz=c*sqrt(t2)
                     t=sz/nc
                     b=bc-(b-sz)
                     d=dc-(d+sz)
                     j=jc-(j+uatan(t,h))
                 endif
         endif
 endif
 return
 end
 !---------------------------------------------------------------------------
 elemental subroutine celbd(mc_qp,elb,eld) !!SJT
 !subroutine celbd(mc,elb,eld) !!SJT: OG
 !
 ! Simultaneous computation of associate complete elliptic integrals
 ! of second kind, B(m) and D(m)
 !
 ! Reference: Fukushima, T (2011) Math. Comp., 80, 1725-1743
 !   Precise and fast computation of general complete elliptic integral
 !   of second kind
 !
 ! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
 !
 real(dp) mc,elk !!SJT
 real(dp),intent(out) :: elb,eld !!SJT
 real(dp) m,mx,kkc,nome
 real(qp),intent(in) :: mc_qp !!SJT
 !
 !real(dp) PIQ,PIHALF!,PI,PIINV !!SJT: removed unused parameters
 real(dp),parameter :: PIQ=0.78539816339744830961566084581988e0_dp
 real(dp),parameter :: PIHALF=1.5707963267948966192313216916398e0_dp
 !parameter (PI=3.1415926535897932384626433832795e0_dp)
 !parameter (PIINV=0.31830988618379067153776752674503e0_dp)
 real(dp) mcold,elbold,eldold
 !save mcold,elbold,eldold !!SJT: disable for thread safety -- variables are initialized at start of each celbd call
 !real(dp) Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16
 real(dp),parameter :: Q1=1.e0_dp/16.e0_dp,Q2=1.e0_dp/32.e0_dp,Q3=21.e0_dp/1024.e0_dp
 real(dp),parameter :: Q4=31.e0_dp/2048.e0_dp,Q5=6257.e0_dp/524288.e0_dp
 real(dp),parameter :: Q6=10293.e0_dp/1048576.e0_dp,Q7=279025.e0_dp/33554432.e0_dp
 real(dp),parameter :: Q8=483127.e0_dp/67108864.e0_dp
 real(dp),parameter :: Q9=435506703.e0_dp/68719476736.e0_dp
 real(dp),parameter :: Q10=776957575.e0_dp/137438953472.e0_dp
 real(dp),parameter :: Q11=22417045555.e0_dp/4398046511104.e0_dp
 real(dp),parameter :: Q12=40784671953.e0_dp/8796093022208.e0_dp
 real(dp),parameter :: Q13=9569130097211.e0_dp/2251799813685248.e0_dp
 real(dp),parameter :: Q14=17652604545791.e0_dp/4503599627370496.e0_dp
 real(dp),parameter :: Q15=523910972020563.e0_dp/144115188075855872.e0_dp
 real(dp),parameter :: Q16=976501268709949.e0_dp/288230376151711744.e0_dp
 !real(dp) K1,K2,K3,K4,K5,K6,K7
 real(dp),parameter :: K1=1.e0_dp/4.e0_dp
 real(dp),parameter :: K2=9.e0_dp/64.e0_dp
 real(dp),parameter :: K3=25.e0_dp/256.e0_dp
 real(dp),parameter :: K4=1225.e0_dp/16384.e0_dp
 real(dp),parameter :: K5=3969.e0_dp/65536.e0_dp
 real(dp),parameter :: K6=53361.e0_dp/1048576.e0_dp
 real(dp),parameter :: K7=184041.e0_dp/4194304.e0_dp
 !real(dp) B1,B2,B3,B4,B5,B6,B7,B8
 real(dp),parameter :: B1=1.e0_dp/2.e0_dp
 real(dp),parameter :: B2=1.e0_dp/16.e0_dp
 real(dp),parameter :: B3=3.e0_dp/128.e0_dp
 real(dp),parameter :: B4=25.e0_dp/2048.e0_dp
 real(dp),parameter :: B5=245.e0_dp/32768.e0_dp
 real(dp),parameter :: B6=1323.e0_dp/262144.e0_dp
 real(dp),parameter :: B7=7623.e0_dp/2097152.e0_dp
 real(dp),parameter :: B8=184041.e0_dp/67108864.e0_dp
 !real(dp) D1,D2,D3,D4,D5,D6,D7,D8
 real(dp),parameter :: D1=1.e0_dp/2.e0_dp
 real(dp),parameter :: D2=3.e0_dp/16.e0_dp
 real(dp),parameter :: D3=15.e0_dp/128.e0_dp
 real(dp),parameter :: D4=175.e0_dp/2048.e0_dp
 real(dp),parameter :: D5=2205.e0_dp/32768.e0_dp
 real(dp),parameter :: D6=14553.e0_dp/262144.e0_dp
 real(dp),parameter :: D7=99099.e0_dp/2097152.e0_dp
 real(dp),parameter :: D8=2760615.e0_dp/67108864.e0_dp
 real(dp) logq2,dkkc,dddc,dele,delb,elk1
 !logical first/.TRUE./    !!!SJT: note that logical variable first is always .true. at the start of each routine call
 logical first !!SJT: initialization within declaration not allowed due to thread safety
 first=.true.  !!SJT: initialize outside of declaration for thread safety
 if(first) then
     first=.FALSE.
         mcold=1.e0_dp
         elbold=PIQ
         eldold=PIQ
 endif
 !m=1.d0-mc !!SJT: original
 m=real(1.0_qp-mc_qp,dp); mc=real(mc_qp,dp) !!SJT
 if(abs(mc-mcold).lt.1.11e-16_dp*mc) then
         elb=elbold
         eld=eldold
 elseif(m.lt.1.11e-16_dp) then
     elb=PIQ
         eld=PIQ
 elseif(mc.lt.1.11e-16_dp) then
     elb=1.e0_dp
         eld=0.3862943611198906188344642429164e0_dp-0.5e0_dp*log(mc)
 elseif(mc.lt.0.1e0_dp) then
     nome=mc*(Q1+mc*(Q2+mc*(Q3+mc*(Q4+mc*(Q5+mc*(Q6 &
         +mc*(Q7+mc*(Q8+mc*(Q9+mc*(Q10+mc*(Q11+mc*(Q12 &
         +mc*(Q13+mc*(Q14+mc*(Q15+mc*Q16))))))))))))))) 
     if(mc.lt.0.01e0_dp) then
         dkkc=mc*(K1+mc*(K2+mc*(K3+mc*(K4+mc*(K5+mc*(K6+mc*K7))))))
         dddc=mc*(D1+mc*(D2+mc*(D3+mc*(D4+mc*(D5+mc*(D6+mc*D7))))))
     else
         mx=mc-0.05e0_dp
 !
 ! (K'-1)/(pi/2)
 !
         dkkc=    0.01286425658832983978282698630501405107893e0_dp &
             +mx*(0.26483429894479586582278131697637750604652e0_dp &
             +mx*(0.15647573786069663900214275050014481397750e0_dp &
             +mx*(0.11426146079748350067910196981167739749361e0_dp &
             +mx*(0.09202724415743445309239690377424239940545e0_dp &
             +mx*(0.07843218831801764082998285878311322932444e0_dp &
             +mx*(0.06935260142642158347117402021639363379689e0_dp &
             +mx*(0.06293203529021269706312943517695310879457e0_dp &
             +mx*(0.05821227592779397036582491084172892108196e0_dp &
             +mx*(0.05464909112091564816652510649708377642504e0_dp &
             +mx*(0.05191068843704411873477650167894906357568e0_dp &
             +mx*(0.04978344771840508342564702588639140680363e0_dp &
             +mx*(0.04812375496807025605361215168677991360500e0_dp &
             ))))))))))))
 !
 ! (K'-E')/(pi/2)
 !
         dddc=    0.02548395442966088473597712420249483947953e0_dp &
             +mx*(0.51967384324140471318255255900132590084179e0_dp &
             +mx*(0.20644951110163173131719312525729037023377e0_dp &
             +mx*(0.13610952125712137420240739057403788152260e0_dp &
             +mx*(0.10458014040566978574883406877392984277718e0_dp &
             +mx*(0.08674612915759188982465635633597382093113e0_dp &
             +mx*(0.07536380269617058326770965489534014190391e0_dp &
             +mx*(0.06754544594618781950496091910264174396541e0_dp &
             +mx*(0.06190939688096410201497509102047998554900e0_dp &
             +mx*(0.05771071515451786553160533778648705873199e0_dp &
             +mx*(0.05451217098672207169493767625617704078257e0_dp &
             +mx*(0.05204028407582600387265992107877094920787e0_dp &
             +mx*(0.05011532514520838441892567405879742720039e0_dp &
             ))))))))))))
     endif
         kkc=1.e0_dp+dkkc
         logq2=-0.5e0_dp*log(nome)
         elk=kkc*logq2
         dele=-dkkc/kkc+logq2*dddc
         elk1=elk-1.e0_dp
         delb=(dele-mc*elk1)/m
         elb=1.e0_dp+delb
         eld=elk1-delb
 elseif(m.le.0.01e0_dp) then
         elb=PIHALF*(B1+m*(B2+m*(B3+m*(B4+m*(B5+m*(B6+m*(B7+m*B8)))))))
         eld=PIHALF*(D1+m*(D2+m*(D3+m*(D4+m*(D5+m*(D6+m*(D7+m*D8)))))))
 elseif(m.le.0.1e0_dp) then
         mx=0.95e0_dp-mc
         elb=     0.790401413584395132310045630540381158921005e0_dp &
             +mx*(0.102006266220019154892513446364386528537788e0_dp &
             +mx*(0.039878395558551460860377468871167215878458e0_dp &
             +mx*(0.021737136375982167333478696987134316809322e0_dp &
             +mx*(0.013960979767622057852185340153691548520857e0_dp &
             +mx*(0.009892518822669142478846083436285145400444e0_dp &
             +mx*(0.007484612400663335676130416571517444936951e0_dp &
             +mx*(0.005934625664295473695080715589652011420808e0_dp &
             +mx*(0.004874249053581664096949448689997843978535e0_dp &
             +mx*(0.004114606930310886136960940893002069423559e0_dp &
             +mx*(0.003550452989196176932747744728766021440856e0_dp &
             +mx*(0.003119229959988474753291950759202798352266e0_dp &
             )))))))))))
         eld=     0.800602040206397047799296975176819811774784e0_dp &
             +mx*(0.313994477771767756849615832867393028789057e0_dp &
             +mx*(0.205913118705551954501930953451976374435088e0_dp &
             +mx*(0.157744346538923994475225014971416837073598e0_dp &
             +mx*(0.130595077319933091909091103101366509387938e0_dp &
             +mx*(0.113308474489758568672985167742047066367053e0_dp &
             +mx*(0.101454199173630195376251916342483192174927e0_dp &
             +mx*(0.0929187842072974367037702927967784464949434e0_dp &
             +mx*(0.0865653801481680871714054745336652101162894e0_dp &
             +mx*(0.0817279846651030135350056216958053404884715e0_dp &
             +mx*(0.0779906657291070378163237851392095284454654e0_dp &
             +mx*(0.075080426851268007156477347905308063808848e0_dp &
             )))))))))))
 elseif(m.le.0.2e0_dp) then
         mx=0.85e0_dp-mc
         elb=     0.80102406445284489393880821604009991524037e0_dp &
             +mx*(0.11069534452963401497502459778015097487115e0_dp &
             +mx*(0.047348746716993717753569559936346358937777e0_dp &
             +mx*(0.028484367255041422845322166419447281776162e0_dp &
             +mx*(0.020277811444003597057721308432225505126013e0_dp &
             +mx*(0.015965005853099119442287313909177068173564e0_dp &
             +mx*(0.013441320273553634762716845175446390822633e0_dp &
             +mx*(0.011871565736951439501853534319081030547931e0_dp &
             +mx*(0.010868363672485520630005005782151743785248e0_dp &
             +mx*(0.010231587232710564565903812652581252337697e0_dp &
             +mx*(0.009849585546666211201566987057592610884309e0_dp &
             +mx*(0.009656606347153765129943681090056980586986e0_dp &
             )))))))))))
         eld=     0.834232667811735098431315595374145207701720e0_dp &
             +mx*(0.360495281619098275577215529302260739976126e0_dp &
             +mx*(0.262379664114505869328637749459234348287432e0_dp &
             +mx*(0.223723944518094276386520735054801578584350e0_dp &
             +mx*(0.206447811775681052682922746753795148394463e0_dp &
             +mx*(0.199809440876486856438050774316751253389944e0_dp &
             +mx*(0.199667451603795274869211409350873244844882e0_dp &
             +mx*(0.204157558868236842039815028663379643303565e0_dp &
             +mx*(0.212387467960572375038025392458549025660994e0_dp &
             +mx*(0.223948914061499360356873401571821627069173e0_dp &
             +mx*(0.238708097425597860161720875806632864507536e0_dp &
             +mx*(0.256707203545463755643710021815937785120030e0_dp &
             )))))))))))
 elseif(m.le.0.3e0_dp) then
         mx=0.75e0_dp-mc
         elb=     0.81259777291992049322557009456643357559904e0_dp &
             +mx*(0.12110961794551011284012693733241967660542e0_dp &
             +mx*(0.057293376831239877456538980381277010644332e0_dp &
             +mx*(0.038509451602167328057004166642521093142114e0_dp &
             +mx*(0.030783430301775232744816612250699163538318e0_dp &
             +mx*(0.027290564934732526869467118496664914274956e0_dp &
             +mx*(0.025916369289445198731886546557337255438083e0_dp &
             +mx*(0.025847203343361799141092472018796130324244e0_dp &
             +mx*(0.026740923539348854616932735567182946385269e0_dp &
             +mx*(0.028464314554825704963640157657034405579849e0_dp &
             +mx*(0.030995446237278954096189768338119395563447e0_dp &
             +mx*(0.034384369179940975864103666824736551261799e0_dp &
             +mx*(0.038738002072493935952384233588242422046537e0_dp &
             ))))))))))))
         eld=     0.873152581892675549645633563232643413901757e0_dp &
             +mx*(0.420622230667770215976919792378536040460605e0_dp &
             +mx*(0.344231061559450379368201151870166692934830e0_dp &
             +mx*(0.331133021818721761888662390999045979071436e0_dp &
             +mx*(0.345277285052808411877017306497954757532251e0_dp &
             +mx*(0.377945322150393391759797943135325823338761e0_dp &
             +mx*(0.427378012464553880508348757311170776829930e0_dp &
             +mx*(0.494671744307822405584118022550673740404732e0_dp &
             +mx*(0.582685115665646200824237214098764913658889e0_dp &
             +mx*(0.695799207728083164790111837174250683834359e0_dp &
             +mx*(0.840018401472533403272555302636558338772258e0_dp &
             +mx*(1.023268503573606060588689738498395211300483e0_dp &
             +mx*(1.255859085136282496149035687741403985044122e0_dp &
             ))))))))))))
 elseif(m.le.0.4e0_dp) then
         mx=0.65e0_dp-mc
         elb=     0.8253235579835158949845697805395190063745e0_dp &
             +mx*(0.1338621160836877898575391383950840569989e0_dp &
             +mx*(0.0710112935979886745743770664203746758134e0_dp &
             +mx*(0.0541784774173873762208472152701393154906e0_dp &
             +mx*(0.0494517449481029932714386586401273353617e0_dp &
             +mx*(0.0502221962241074764652127892365024315554e0_dp &
             +mx*(0.0547429131718303528104722303305931350375e0_dp &
             +mx*(0.0627462579270016992000788492778894700075e0_dp &
             +mx*(0.0746698810434768864678760362745179321956e0_dp &
             +mx*(0.0914808451777334717996463421986810092918e0_dp &
             +mx*(0.1147050921109978235104185800057554574708e0_dp &
             +mx*(0.1465711325814398757043492181099197917984e0_dp &
             +mx*(0.1902571373338462844225085057953823854177e0_dp &
             ))))))))))))
         eld=     0.9190270392420973478848471774160778462738e0_dp &
             +mx*(0.5010021592882475139767453081737767171354e0_dp &
             +mx*(0.4688312705664568629356644841691659415972e0_dp &
             +mx*(0.5177142277764000147059587510833317474467e0_dp &
             +mx*(0.6208433913173031070711926900889045286988e0_dp &
             +mx*(0.7823643937868697229213240489900179142670e0_dp &
             +mx*(1.0191145350761029126165253557593691585239e0_dp &
             +mx*(1.3593452027484960522212885423056424704073e0_dp &
             +mx*(1.8457173023588279422916645725184952058635e0_dp &
             +mx*(2.5410717031539207287662105618152273788399e0_dp &
             +mx*(3.5374046552080413366422791595672470037341e0_dp &
             +mx*(4.9692960029774259303491034652093672488707e0_dp &
             +mx*(7.0338228700300311264031522795337352226926e0_dp &
             +mx*(10.020043225034471401553194050933390974016e0_dp &
             )))))))))))))
 elseif(m.le.0.5e0_dp) then
         mx=0.55e0_dp-mc
         elb=     0.8394795702706129706783934654948360410325e0_dp &
             +mx*(0.1499164403063963359478614453083470750543e0_dp &
             +mx*(0.0908319358194288345999005586556105610069e0_dp &
             +mx*(0.0803470334833417864262134081954987019902e0_dp &
             +mx*(0.0856384405004704542717663971835424473169e0_dp &
             +mx*(0.1019547259329903716766105911448528069506e0_dp &
             +mx*(0.1305748115336160150072309911623351523284e0_dp &
             +mx*(0.1761050763588499277679704537732929242811e0_dp &
             +mx*(0.2468351644029554468698889593583314853486e0_dp &
             +mx*(0.3564244768677188553323196975301769697977e0_dp &
             +mx*(0.5270025622301027434418321205779314762241e0_dp &
             +mx*(0.7943896342593047502260866957039427731776e0_dp &
             +mx*(1.2167625324297180206378753787253096783993e0_dp &
             ))))))))))))
         eld=     0.9744043665463696730314687662723484085813e0_dp &
             +mx*(0.6132468053941609101234053415051402349752e0_dp &
             +mx*(0.6710966695021669963502789954058993004082e0_dp &
             +mx*(0.8707276201850861403618528872292437242726e0_dp &
             +mx*(1.2295422312026907609906452348037196571302e0_dp &
             +mx*(1.8266059675444205694817638548699906990301e0_dp &
             +mx*(2.8069345309977627400322167438821024032409e0_dp &
             +mx*(4.4187893290840281339827573139793805587268e0_dp &
             +mx*(7.0832360574787653249799018590860687062869e0_dp &
             +mx*(11.515088120557582942290563338274745712174e0_dp &
             +mx*(18.931511185999274639516732819605594455165e0_dp &
             +mx*(31.411996938204963878089048091424028309798e0_dp &
             +mx*(52.520729454575828537934780076768577185134e0_dp &
             +mx*(88.384854735065298062125622417251073520996e0_dp &
             +mx*(149.56637449398047835236703116483092644714e0_dp &
             +mx*(254.31790843104117434615624121937495622372e0_dp &
             )))))))))))))))
 elseif(m.le.0.6e0_dp) then
         mx=0.45e0_dp-mc
         elb=     0.8554696151564199914087224774321783838373e0_dp &
             +mx*(0.1708960726897395844132234165994754905373e0_dp &
             +mx*(0.1213352290269482260207667564010437464156e0_dp &
             +mx*(0.1282018835749474096272901529341076494573e0_dp &
             +mx*(0.1646872814515275597348427294090563472179e0_dp &
             +mx*(0.2374189087493817423375114793658754489958e0_dp &
             +mx*(0.3692081047164954516884561039890508294508e0_dp &
             +mx*(0.6056587338479277173311618664015401963868e0_dp &
             +mx*(1.0337055615578127436826717513776452476106e0_dp &
             +mx*(1.8189884893632678849599091011718520567105e0_dp &
             +mx*(3.2793776512738509375806561547016925831128e0_dp &
             +mx*(6.0298883807175363312261449542978750456611e0_dp &
             +mx*(11.269796855577941715109155203721740735793e0_dp &
             +mx*(21.354577850382834496786315532111529462693e0_dp &
             )))))))))))))
         eld=     1.04345529511513353426326823569160142342838e0_dp &
             +mx*(0.77962572192850485048535711388072271372632e0_dp &
             +mx*(1.02974236093206758187389128668777397528702e0_dp &
             +mx*(1.62203722341135313022433907993860147395972e0_dp &
             +mx*(2.78798953118534762046989770119382209443756e0_dp &
             +mx*(5.04838148737206914685643655935236541332892e0_dp &
             +mx*(9.46327761194348429539987572314952029503864e0_dp &
             +mx*(18.1814899494276679043749394081463811247757e0_dp &
             +mx*(35.5809805911791687037085198750213045708148e0_dp &
             +mx*(70.6339354619144501276254906239838074917358e0_dp &
             +mx*(141.828580083433059297030133195739832297859e0_dp &
             +mx*(287.448751250132166257642182637978103762677e0_dp &
             +mx*(587.115384649923076181773192202238389711345e0_dp &
             +mx*(1207.06543522548061603657141890778290399603e0_dp &
             +mx*(2495.58872724866422273012188618178997342537e0_dp &
             +mx*(5184.69242939480644062471334944523925163600e0_dp &
             +mx*(10817.2133369041327524988910635205356016939e0_dp &
             ))))))))))))))))
 elseif(m.le.0.7e0_dp) then
         mx=0.35e0_dp-mc
         elb=     0.8739200618486431359820482173294324246058e0_dp &
             +mx*(0.1998140574823769459497418213885348159654e0_dp &
             +mx*(0.1727696158780152128147094051876565603862e0_dp &
             +mx*(0.2281069132842021671319791750725846795701e0_dp &
             +mx*(0.3704681411180712197627619157146806221767e0_dp &
             +mx*(0.6792712528848205545443855883980014994450e0_dp &
             +mx*(1.3480084966817573020596179874311042267679e0_dp &
             +mx*(2.8276709768538207038046918250872679902352e0_dp &
             +mx*(6.1794682501239140840906583219887062092430e0_dp &
             +mx*(13.935686010342811497608625663457407447757e0_dp &
             +mx*(32.218929281059722026322932181420383764028e0_dp &
             +mx*(76.006962959226101026399085304912635262362e0_dp &
             +mx*(182.32144908775406957609058046006949657416e0_dp &
             +mx*(443.51507644112648158679360783118806161062e0_dp &
             +mx*(1091.8547229028388292980623647414961662223e0_dp &
             +mx*(2715.7658664038195881056269799613407111521e0_dp &
             )))))))))))))))
         eld=     1.13367833657573316571671258513452768536080e0_dp &
             +mx*(1.04864317372997039116746991765351150490010e0_dp &
             +mx*(1.75346504119846451588826580872136305225406e0_dp &
             +mx*(3.52318272680338551269021618722443199230946e0_dp &
             +mx*(7.74947641381397458240336356601913534598302e0_dp &
             +mx*(17.9864500558507330560532617743406294626849e0_dp &
             +mx*(43.2559163462326133313977294448984936591235e0_dp &
             +mx*(106.681534454096017031613223924991564429656e0_dp &
             +mx*(268.098486573117433951562111736259672695883e0_dp &
             +mx*(683.624114850289804796762005964155730439745e0_dp &
             +mx*(1763.49708521918740723028849567007874329637e0_dp &
             +mx*(4592.37475383116380899419201719007475759114e0_dp &
             +mx*(12053.4410190488892782190764838488156555734e0_dp &
             +mx*(31846.6630207420816960681624497373078887317e0_dp &
             +mx*(84621.2213590568080177035346867495326879117e0_dp &
             +mx*(225956.423182907889987641304430180593010940e0_dp &
             +mx*(605941.517281758859958050194535269219533685e0_dp &
             +mx*(1.63108259953926832083633544697688841456604e6_dp &
             )))))))))))))))))
 elseif(m.le.0.8e0_dp) then
         mx=0.25e0_dp-mc
         elb=     0.895902820924731621258525533131864225704e0_dp &
             +mx*(0.243140003766786661947749288357729051637e0_dp &
             +mx*(0.273081875594105531575351304277604081620e0_dp &
             +mx*(0.486280007533573323895498576715458103274e0_dp &
             +mx*(1.082747437228230914750752674136983406683e0_dp &
             +mx*(2.743445290986452500459431536349945437824e0_dp &
             +mx*(7.555817828670234627048618342026400847824e0_dp &
             +mx*(22.05194082493752427472777448620986154515e0_dp &
             +mx*(67.15640644740229407624192175802742979626e0_dp &
             +mx*(211.2722537881770961691291434845898538537e0_dp &
             +mx*(681.9037843053270682273212958093073895805e0_dp &
             +mx*(2246.956231592536516768812462150619631201e0_dp &
             +mx*(7531.483865999711792004783423815426725079e0_dp &
             +mx*(25608.51260130241579018675054866136922157e0_dp &
             +mx*(88140.74740089604971425934283371277143256e0_dp &
             +mx*(306564.4242098446591430938434419151070722e0_dp &
             +mx*(1.076036077811072193752770590363885180738e6_dp &
             +mx*(3.807218502573632648224286313875985190526e6_dp &
             +mx*(1.356638224422139551020110323739879481042e7_dp &
             ))))))))))))))))))
         eld=     1.26061282657491161418014946566845780315983e0_dp &
             +mx*(1.54866563808267658056930177790599939977154e0_dp &
             +mx*(3.55366941187160761540650011660758187283401e0_dp &
             +mx*(9.90044467610439875577300608183010716301714e0_dp &
             +mx*(30.3205666174524719862025105898574414438275e0_dp &
             +mx*(98.1802586588830891484913119780870074464833e0_dp &
             +mx*(329.771010434557055036273670551546757245808e0_dp &
             +mx*(1136.65598974289039303581967838947708073239e0_dp &
             +mx*(3993.83433574622979757935610692842933356144e0_dp &
             +mx*(14242.7295865552708506232731633468180669284e0_dp &
             +mx*(51394.7572916887209594591528374806790960057e0_dp &
             +mx*(187246.702914623152141768788258141932569037e0_dp &
             +mx*(687653.092375389902708761221294282367947659e0_dp &
             +mx*(2.54238553565398227033448846432182516906624e6_dp &
             +mx*(9.45378121934749027243313241962076028066811e6_dp &
             +mx*(3.53283630179709170835024033154326126569613e7_dp &
             +mx*(1.32593262383393014923560730485845833322771e8_dp &
             +mx*(4.99544968184054821463279808395426941549833e8_dp &
             +mx*(1.88840934729443872364972817525484292678543e9_dp &
             +mx*(7.16026753447893719179055010636502508063102e9_dp &
             +mx*(2.72233079469633962247554894093591262281929e10_dp &
         ))))))))))))))))))))
 elseif(m.le.0.85e0_dp) then
         mx=0.175e0_dp-mc
         elb=     0.915922052601931494319853880201442948834592e0_dp &
             +mx*(0.294714252429483394379515488141632749820347e0_dp &
             +mx*(0.435776709264636140422971598963772380161131e0_dp &
             +mx*(1.067328246493644238508159085364429570207744e0_dp &
             +mx*(3.327844118563268085074646976514979307993733e0_dp &
             +mx*(11.90406004445092906188837729711173326621810e0_dp &
             +mx*(46.47838820224626393512400481776284680677096e0_dp &
             +mx*(192.7556002578809476962739389101964074608802e0_dp &
             +mx*(835.3356299261900063712302517586717381557137e0_dp &
             +mx*(3743.124548343029102644419963712353854902019e0_dp &
             +mx*(17219.07731004063094108708549153310467326395e0_dp &
             +mx*(80904.60401669850158353080543152212152282878e0_dp &
             +mx*(386808.3292751742460123683674607895217760313e0_dp &
             +mx*(1.876487670110449342170327796786290400635732e6_dp &
             +mx*(9.216559908641567755240142886998737950775910e6_dp &
             ))))))))))))))
         eld=     1.402200569110579095046054435635136986038164e0_dp &
             +mx*(2.322205897861749446534352741005347103992773e0_dp &
             +mx*(7.462158366466719682730245467372788273333992e0_dp &
             +mx*(29.43506890797307903104978364254987042421285e0_dp &
             +mx*(128.1590924337895775262509354898066132182429e0_dp &
             +mx*(591.0807036911982326384997979640812493154316e0_dp &
             +mx*(2830.546229607726377048576057043685514661188e0_dp &
             +mx*(13917.76431889392229954434840686741305556862e0_dp &
             +mx*(69786.10525163921228258055074102587429394212e0_dp &
             +mx*(355234.1420341879634781808899208309503519936e0_dp &
             +mx*(1.830019186413931053503912913904321703777885e6_dp &
             +mx*(9.519610812032515607466102200648641326190483e6_dp &
             +mx*(4.992086875574849453986274042758566713803723e7_dp &
             +mx*(2.635677009826023473846461512029006874800883e8_dp &
             +mx*(1.399645765120061118824228996253541612110338e9_dp &
             +mx*(7.469935792837635004663183580452618726280406e9_dp &
             +mx*(4.004155595835610574316003488168804738481448e10_dp &
             +mx*(2.154630668144966654449602981243932210695662e11_dp &
             )))))))))))))))))
 else
     mx=0.125e0_dp-mc
         elb=     0.931906061029524827613331428871579482766771e0_dp &
             +mx*(0.348448029538453860999386797137074571589376e0_dp &
             +mx*(0.666809178846938247558793864839434184202736e0_dp &
             +mx*(2.210769135708128662563678717558631573758222e0_dp &
             +mx*(9.491765048913406881414290930355300611703187e0_dp &
             +mx*(47.09304791027740853381457907791343619298913e0_dp &
             +mx*(255.9200460211233087050940506395442544885608e0_dp &
             +mx*(1480.029532675805407554800779436693505109703e0_dp &
             +mx*(8954.040904734313578374783155553041065984547e0_dp &
             +mx*(56052.48220982686949967604699243627330816542e0_dp &
             +mx*(360395.7241626000916973524840479780937869149e0_dp &
             +mx*(2.367539415273216077520928806581689330885103e6_dp &
             +mx*(1.582994957277684102454906900025484391190264e7_dp &
             +mx*(1.074158093278511100137056972128875270067228e8_dp &
             +mx*(7.380585460239595691878086073095523043390649e8_dp &
             +mx*(5.126022002555101496684687154904781856830296e9_dp &
             +mx*(3.593534065502416588712409180013118409428367e10_dp &
             +mx*(2.539881257612812212004146637239987308133582e11_dp &
             +mx*(1.808180007145359569674767150594344316702507e12_dp &
         ))))))))))))))))))
         eld=     1.541690112721819084362258323861459983048179e0_dp &
             +mx*(3.379176214579645449453938918349243359477706e0_dp &
             +mx*(14.94058385670236671625328259137998668324435e0_dp &
             +mx*(81.91773929235074880784578753539752529822986e0_dp &
             +mx*(497.4900546551479866036061853049402721939835e0_dp &
             +mx*(3205.184010234846235275447901572262470252768e0_dp &
             +mx*(21457.32237355321925571253220641357074594515e0_dp &
             +mx*(147557.0156564174712105689758692929775004292e0_dp &
             +mx*(1.035045290185256525452269053775273002725343e6_dp &
             +mx*(7.371922334832212125197513363695905834126154e6_dp &
             +mx*(5.314344395142401141792228169170505958906345e7_dp &
             +mx*(3.868823475795976312985118115567305767603128e8_dp &
             +mx*(2.839458401528033778425531336599562337200510e9_dp &
             +mx*(2.098266122943898941547136470383199468548861e10_dp &
             +mx*(1.559617754017662417944194874282275405676282e11_dp &
             +mx*(1.165096220419884791236699872205721392201682e12_dp &
             +mx*(8.742012983013913804987431275193291316808766e12_dp &
             +mx*(6.584725462672366918676967847406180155459650e13_dp &
             +mx*(4.976798737062434393396993620379481464465749e14_dp &
             +mx*(3.773018634056605404718444239040628892506293e15_dp &
             +mx*(2.868263194837819660109735981973458220407767e16_dp &
         ))))))))))))))))))))
 endif
 mcold=mc
 elbold=elb
 eldold=eld
 return
 end
 !---------------------------------------------------------------------------
 elemental subroutine celbdj(nc0,mc0_qp,celb,celd,celj)
 !subroutine celbdj(nc0,mc0,celb,celd,celj) !!SJT: OG
 !
 ! Simultaneous computation of associate complete elliptic integrals
 ! of third kind, B(m), D(m), and J(n|m)
 !
 ! Reference: Fukushima, T, (2013) J. Comp. Appl. Math., 253, 142-157
 !     Fast computation of a general complete elliptic integral
 !     of third kind by half and double argument transformations
 !
 ! Restricted Version with Assumptions: 0 < mc, 0 < nc
 !
 ! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
 !
 real(dp),intent(in) :: nc0
 real(dp),intent(out) :: celb,celd,celj !!SJT
 logical flag
 !integer(isp) IMAX,i,is,ie !!SJT: OG
 integer(isp) i,is,ie !!SJT
 integer(isp),parameter :: IMAX=40_isp
 real(dp) y(0_isp:IMAX),x(0_isp:IMAX),c(0_isp:IMAX),d(0_isp:IMAX),a(0_isp:IMAX)
 real(dp) mc,mc0,nc,m,n,celk,yi,ye,dj,m1,kc0,temp !!SJT
 real(qp),intent(in) :: mc0_qp !!SJT
 real(qp) :: mc_qp !!SJT
 !real(dp) PIHALF,EPS,THIRD
 real(dp),parameter :: PIHALF=1.5707963267948966e0_dp
 real(dp),parameter :: EPS=1.11e-16_dp
 real(dp),parameter :: THIRD=1.e0_dp/3.e0_dp
 real(dp) B(IMAX)
 !logical first /.TRUE./  !!!SJT: note that logical variable first is .true. at the start of each routine call
 !save B,first !!SJT: disable for thread safety
 logical first !!SJT: initialization cannot be in declaration for thread-safe calculations
 first=.true. !!SJT: initialization outside of declaration
 mc0=real(mc0_qp,dp) !!SJT
 !if(mc0.le.0.d0) then !!SJT: does not occur in practice
 !    write(*,*) "(celbdj) Out of domain: mc <= 0"
 !    return
 !endif
 !if(nc0.le.0.d0) then
 !    write(*,*) "(celbdj) Out of domain: nc <= 0"
 !    return
 !endif
 if(mc0.lt.1.e0_dp) then
     !mc=mc0;nc=nc0
     mc_qp=mc0_qp; nc=nc0 !!SJT
 elseif(mc0.gt.1.e0_dp) then
     !mc=1.e0_dp/mc0;nc=nc0*mc !!SJT: OG
     mc_qp=1.0_qp/mc0_qp; nc=nc0*real(mc_qp,dp) !!SJT
 else
     celb=PIHALF*0.5e0_dp
     celd=PIHALF*0.5e0_dp
     celj=PIHALF/(nc0+sqrt(nc0))
     return
 endif
 if(first) then
     first=.FALSE.
     do i=1_isp,IMAX
         B(i)=1.e0_dp/dble(2_isp*i+1_isp)
     enddo
 endif
 mc=real(mc_qp,dp) !!SJT
 !call celbd(mc,celb,celd) !!SJT: OG
 call celbd(mc_qp,celb,celd) !!SJT
 if(nc.eq.mc) then
     celj=celb/mc
     goto 4  
 endif
 !m=1.d0-mc; n=1.d0-nc !!SJT: original
 m=real(1.0_qp-mc_qp,dp); n=1.e0_dp-nc !!SJT
 flag=nc.lt.mc.or.(n*nc).gt.(nc-mc)
 if(flag) then
     y(0_isp)=(nc-mc)/(nc*m)
 else
     y(0_isp)=n/m
 endif
 is=0_isp
 if(y(0_isp).gt.0.5d0) then
     x(0_isp)=1.e0_dp-y(0_isp)
     do i=0_isp,IMAX-1_isp !!SJT
     !do i=0,IMAX !!SJT: OG
         c(i)=sqrt(x(i))
         d(i)=sqrt(mc+m*x(i))
         x(i+1_isp)=(c(i)+d(i))/(1.e0_dp+d(i))
         if(x(i+1_isp).gt.0.5e0_dp) then
             y(i+1_isp)=y(i)/((1.e0_dp+c(i))*(1.e0_dp+d(i)))
             is=i+1_isp
             goto 1
         endif
         y(i+1_isp)=1.e0_dp-x(i+1_isp)
     enddo
     !write(*,"(a30,i5)") "(celbdj) No Conv. x-Transf. i=",i !!SJT: does not occur in practice
     return
 endif
 1 continue
 do i=is,IMAX-1_isp !!SJT
 !do i=is,IMAX !!SJT: OG
     c(i)=sqrt(1.e0_dp-y(i))
     d(i)=sqrt(1.e0_dp-m*y(i))
     y(i+1_isp)=y(i)/((1.e0_dp+c(i))*(1.e0_dp+d(i)))
     if(abs(y(i+1_isp)).lt.0.325e0_dp) goto 2
 enddo
 !write(*,"(a30,i5)") "(celbdj) No Conv. y-Transf. i=",i !!SJT: does not occur in practice
 return
 2 continue
 ie=i+1_isp
 ye=y(ie)
 celk=celb+celd
 a(0_isp)=celd
 celj=a(0_isp)
 yi=ye
 a(1_isp)=((1.e0_dp+2.e0_dp*m)*celd-celb)*THIRD
 dj=a(1_isp)*yi
 i=1_isp
 if(abs(dj).lt.EPS*abs(celj)) goto 3
 celj=celj+dj
 m1=1.e0_dp+m
 do i=2_isp,IMAX
     yi=yi*ye
     a(i)=(1.e0_dp-B(i))*m1*a(i-1_isp)-(1.e0_dp-2.e0_dp*B(i))*m*a(i-2_isp)
     dj=a(i)*yi
     if(abs(dj).lt.EPS*abs(celj)) goto 3
     celj=celj+dj
 enddo
 !write(*,"(a30,i5)") "(celbdj) No Conv. Series. i=",i !!SJT: does not occur in practice
 return
 3 continue
 do i=ie-1_isp,0_isp,-1_isp
     celj=(2.e0_dp*(c(i)+d(i))*celj-y(i)*celk)/(c(i)*d(i)*(1.e0_dp+c(i))*(1.e0_dp+d(i)))
 enddo
 if(flag) then
     celj=(nc*celk-mc*celj)/(nc*nc)
 endif
 4 continue
 if(mc0.gt.1.e0_dp) then
     kc0=sqrt(mc0)
     temp=celb
     celb=celd/kc0
     celd=temp/kc0
     celj=celj/(mc0*kc0)
 endif
 !write(*,"(a30,1p3e15.7)") "(celbdj) nc0,mc0,celj=",nc0,mc0,celj
 return
 end
 !---------------------------------------------------------------------------
 elemental subroutine elsbdj(s0,n,mc_qp,b,d,j) !!SJT
 !subroutine elsbdj(s0,n,mc,b,d,j) !!SJT: OG
 !
 ! Simultaneous computation of associate elliptic integrals
 ! of third kind, B(phi|m), D(phi|m), and J(phi,n|m)
 ! by using the half/double argument transformation of sn functions
 !
 ! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !   Precise and fast computation of a general incomplete elliptic integral
 !   of third kind by half and double argument transformations
 !
 real(dp),intent(in) :: s0,n !!SJT
 real(dp),intent(out) :: b,d,j !!SJT
 !real(dp) m,mc,h,del,s,y,c,sy,t !!SJT
 real(dp) m,h,del,s,y,c,sy,t !!SJT
 real(qp),intent(in) :: mc_qp !!SJT
 real(dp) yy(11_isp),ss(11_isp),cd(11_isp)
 !real(dp) serj,uatan !!SJT: due to module definition, this routine already has access to serj and uatan functions
 integer(isp) i,k
 
 ! write(*,*) "(elsj) s0,n,mc=",s0,n,mc
 
 !m=1.d0-mc !!SJT: original
 m=real(1.0_qp-mc_qp,dp); !mc=real(mc_qp,dp) !!SJT
 !kc=mc**0.5e0_dp; m=(1.e0_dp+kc)*(1.e0_dp-kc) !!SJT: to reduce truncation error when mc is close to unity
 h=n*(1.e0_dp-n)*(n-m)
 !del=0.04094e0_dp-0.00652e0_dp*m
 !del=0.03429e0_dp
 !del=0.02448e0_dp  ! JA
 del=0.01622e0_dp   ! J9
 !del=0.00969e0_dp  ! J8
 !del=0.00500e0_dp  ! J7
 !del=2.073e-3_dp    ! J6
 !del=6.037e-4_dp    ! J5
 
 s=s0
 y=s*s
 if(y.lt.del) then
         call serbd(y,m,b,d)
         b=s*b
         d=s*y*d
         j=s*serj(y,n,m)
 ! write(*,"(a20,1p3e10.2)") "(elsbdj) b,d,j=",b,d,j
         return
 endif
 yy(1_isp)=y
 ss(1_isp)=s
 ! write(*,"(a20,i10,1pe10.2)") "(elsbdj) i,y=",1,y
 do i=1_isp,10_isp
     c=sqrt(1.e0_dp-y)
     d=sqrt(1.e0_dp-m*y)
     y=y/((1.e0_dp+c)*(1.e0_dp+d))
         yy(i+1_isp)=y
         ss(i+1_isp)=sqrt(y)
         cd(i)=c*d
 ! write(*,"(a20,i10,1pe10.2)") "(elsbdj) i,y=",i+1,y
         if(y.lt.del) then
                 goto 1
         endif
 enddo
 !write(*,*) "(elsbdj) too many iterations: s0,m=",s0,m !SJT: does not occur in practice
 1 continue
 ! write(*,"(a20,i10,1pe10.2)") "(elsbdj) i,y=",i+1,y
 call serbd(y,m,b,d)
 b=ss(i+1_isp)*b
 d=ss(i+1_isp)*y*d
 j=ss(i+1_isp)*serj(y,n,m)
 do k=i,1_isp,-1_isp
         sy=ss(k)*yy(k+1_isp)
         t=sy/(1.e0_dp-n*(yy(k)-yy(k+1_isp)*cd(k)))
         b=2.e0_dp*b-sy
         d=d+(d+sy)
         j=j+(j+uatan(t,h))
 enddo
 ! write(*,"(a20,1p3e10.2)") "(elsbdj) b,d,j=",b,d,j
 return
 end
 !---------------------------------------------------------------------------
 elemental subroutine serbd(y,m,b,d)
 !
 ! Simultaneous computation of associate elliptic integrals,
 ! B(phi|m) and D(phi|m), for small arguments by the series expansion
 !
 ! Reference: Fukushima, T (2012) J. Comp. Appl. Math., 235, 4140-4148
 !   Precise and fast computation of general incomplete elliptic integral
 !   of second kind by half and double argument transformations
 !
 real(dp),intent(in) :: y,m !!SJT
 real(dp),intent(out) :: b,d !!SJT
 !
 !real(dp) F1,F2,F3,F4 !!SJT: OG
 real(dp) F1,F2,F3,F4,F5,F6,F7,F8,F9,FA,FB !!SJT: non-parameter variables
 !real(dp) F10,F20,F21,F30,F31,F40,F41,F42
 !real(dp) F5,F50,F51,F52,F6,F60,F61,F62,F63
 !real(dp) F7,F70,F71,F72,F73,F8,F80,F81,F82,F83,F84
 !real(dp) F9,F90,F91,F92,F93,F94
 !real(dp) FA,FA0,FA1,FA2,FA3,FA4,FA5
 !real(dp) FB,FB0,FB1,FB2,FB3,FB4,FB5
 real(dp),parameter :: F10=1.e0_dp/6.e0_dp
 real(dp),parameter :: F20=3.e0_dp/40.e0_dp
 real(dp),parameter :: F21=2.e0_dp/40.e0_dp
 real(dp),parameter :: F30=5.e0_dp/112.e0_dp
 real(dp),parameter :: F31=3.e0_dp/112.e0_dp
 real(dp),parameter :: F40=35.e0_dp/1152.e0_dp
 real(dp),parameter :: F41=20.e0_dp/1152.e0_dp
 real(dp),parameter :: F42=18.e0_dp/1152.e0_dp
 real(dp),parameter :: F50=63.e0_dp/2816.e0_dp
 real(dp),parameter :: F51=35.e0_dp/2816.e0_dp
 real(dp),parameter :: F52=30.e0_dp/2816.e0_dp
 real(dp),parameter :: F60=231.e0_dp/13312.e0_dp
 real(dp),parameter :: F61=126.e0_dp/13312.e0_dp
 real(dp),parameter :: F62=105.e0_dp/13312.e0_dp
 real(dp),parameter :: F63=100.e0_dp/13312.e0_dp
 real(dp),parameter :: F70=429.e0_dp/30720.e0_dp
 real(dp),parameter :: F71=231.e0_dp/30720.e0_dp
 real(dp),parameter :: F72=189.e0_dp/30720.e0_dp
 real(dp),parameter :: F73=175.e0_dp/30720.e0_dp
 real(dp),parameter :: F80=6435.e0_dp/557056.e0_dp
 real(dp),parameter :: F81=3432.e0_dp/557056.e0_dp
 real(dp),parameter :: F82=2722.e0_dp/557056.e0_dp
 real(dp),parameter :: F83=2520.e0_dp/557056.e0_dp
 real(dp),parameter :: F84=2450.e0_dp/557056.e0_dp
 real(dp),parameter :: F90=12155.e0_dp/1245184.e0_dp
 real(dp),parameter :: F91=6435.e0_dp/1245184.e0_dp
 real(dp),parameter :: F92=5148.e0_dp/1245184.e0_dp
 real(dp),parameter :: F93=4620.e0_dp/1245184.e0_dp
 real(dp),parameter :: F94=4410.e0_dp/1245184.e0_dp
 real(dp),parameter :: FA0=46189.e0_dp/5505024.e0_dp
 real(dp),parameter :: FA1=24310.e0_dp/5505024.e0_dp
 real(dp),parameter :: FA2=19305.e0_dp/5505024.e0_dp
 real(dp),parameter :: FA3=17160.e0_dp/5505024.e0_dp
 real(dp),parameter :: FA4=16170.e0_dp/5505024.e0_dp
 real(dp),parameter :: FA5=15876.e0_dp/5505024.e0_dp
 real(dp),parameter :: FB0=88179.e0_dp/12058624.e0_dp
 real(dp),parameter :: FB1=46189.e0_dp/12058624.e0_dp
 real(dp),parameter :: FB2=36465.e0_dp/12058624.e0_dp
 real(dp),parameter :: FB3=32175.e0_dp/12058624.e0_dp
 real(dp),parameter :: FB4=30030.e0_dp/12058624.e0_dp
 real(dp),parameter :: FB5=29106.e0_dp/12058624.e0_dp
 !real(dp) A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB
 real(dp),parameter :: A1=3.e0_dp/5.e0_dp
 real(dp),parameter :: A2=5.e0_dp/7.e0_dp
 real(dp),parameter :: A3=7.e0_dp/9.e0_dp
 real(dp),parameter :: A4=9.e0_dp/11.e0_dp
 real(dp),parameter :: A5=11.e0_dp/13.e0_dp
 real(dp),parameter :: A6=13.e0_dp/15.e0_dp
 real(dp),parameter :: A7=15.e0_dp/17.e0_dp
 real(dp),parameter :: A8=17.e0_dp/19.e0_dp
 real(dp),parameter :: A9=19.e0_dp/21.e0_dp
 real(dp),parameter :: AA=21.e0_dp/23.e0_dp
 real(dp),parameter :: AB=23.e0_dp/25.e0_dp
 real(dp) B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB
 !real(dp) D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,DA,DB !!SJT: OG
 real(dp) D1,D2,D3,D4,D5,D6,D7,D8,D9,DA,DB !!SJT: non-parameter variables
 real(dp),parameter :: D0=1.e0_dp/3.e0_dp
 !
 !        write(*,*) "(serbd) y,m=",y,m
 F1=F10+m*F10
 F2=F20+m*(F21+m*F20)
 F3=F30+m*(F31+m*(F31+m*F30))
 F4=F40+m*(F41+m*(F42+m*(F41+m*F40)))
 F5=F50+m*(F51+m*(F52+m*(F52+m*(F51+m*F50))))
 F6=F60+m*(F61+m*(F62+m*(F63+m*(F62+m*(F61+m*F60)))))
 F7=F70+m*(F71+m*(F72+m*(F73+m*(F73+m*(F72+m*(F71+m*F70))))))
 F8=F80+m*(F81+m*(F82+m*(F83+m*(F84+m*(F83+m*(F82+m*(F81+m*F80)))))))
 F9=F90+m*(F91+m*(F92+m*(F93+m*(F94+m*(F94+m*(F93 &
     +m*(F92+m*(F91+m*F90))))))))
 FA=FA0+m*(FA1+m*(FA2+m*(FA3+m*(FA4+m*(FA5+m*(FA4 &
     +m*(FA3+m*(FA2+m*(FA1+m*FA0)))))))))
 FB=FB0+m*(FB1+m*(FB2+m*(FB3+m*(FB4+m*(FB5+m*(FB5+ &
     m*(FB4+m*(FB3+m*(FB2+m*(FB1+m*FB0))))))))))
 !
 D1=F1*A1
 D2=F2*A2
 D3=F3*A3
 D4=F4*A4
 D5=F5*A5
 D6=F6*A6
 D7=F7*A7
 D8=F8*A8
 D9=F9*A9
 DA=FA*AA
 DB=FB*AB
 d=D0+y*(D1+y*(D2+y*(D3+y*(D4+y*(D5+y*(D6+y*(D7+y*(D8 &
     +y*(D9+y*(DA+y*DB))))))))))
 B1=F1-D0
 B2=F2-D1
 B3=F3-D2
 B4=F4-D3
 B5=F5-D4
 B6=F6-D5
 B7=F7-D6
 B8=F8-D7
 B9=F9-D8
 BA=FA-D9
 BB=FB-DA
 b=1.e0_dp+y*(B1+y*(B2+y*(B3+y*(B4+y*(B5+y*(B6+y*(B7+y*(B8 &
     +y*(B9+y*(BA+y*BB))))))))))
 return
 end
 !---------------------------------------------------------------------------
 elemental real(dp) function serj(y,n,m)
 !
 ! Computation of associate elliptic integral J(phi,n|m)
 ! for small arguments by the series expansion
 !
 ! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !   Precise and fast computation of a general incomplete elliptic integral
 !   of third kind by half and double argument transformations
 !
 real(dp),intent(in) :: y,n,m !!SJT
 
 real(dp) J1,J2,J3,J4,J5,J6,J7,J8,J9,JA
 
 !real(dp) J100,J200,J201,J210,J300,J301,J302,J310,J311,J320
 !real(dp) J400,J401,J402,J403,J410,J411,J412,J420,J421,J430
 !real(dp) J500,J501,J502,J503,J504,J510,J511,J512,J513,J520
 !real(dp) J521,J522,J530,J531,J540
 !real(dp) J600,J601,J602,J603,J604,J605,J610,J611,J612,J613,J614
 !real(dp) J620,J621,J622,J623,J630,J631,J632,J640,J641,J650
 !real(dp) J700,J701,J702,J703,J704,J705,J706
 !real(dp) J710,J711,J712,J713,J714,J715,J720,J721,J722,J723,J724
 !real(dp) J730,J731,J732,J733,J740,J741,J742,J750,J751,J760
 !real(dp) J800,J801,J802,J803,J804,J805,J806,J807
 !real(dp) J810,J811,J812,J813,J814,J815,J816
 !real(dp) J820,J821,J822,J823,J824,J825,J830,J831,J832,J833,J834
 !real(dp) J840,J841,J842,J843,J850,J851,J852,J860,J861,J870
 !real(dp) J900,J901,J902,J903,J904,J905,J906,J907,J908
 !real(dp) J910,J911,J912,J913,J914,J915,J916,J917
 !real(dp) J920,J921,J922,J923,J924,J925,J926
 !real(dp) J930,J931,J932,J933,J934,J935,J940,J941,J942,J943,J944
 !real(dp) J950,J951,J952,J953,J960,J961,J962,J970,J971,J980
 !real(dp) JA00,JA01,JA02,JA03,JA04,JA05,JA06,JA07,JA08,JA09
 !real(dp) JA10,JA11,JA12,JA13,JA14,JA15,JA16,JA17,JA18
 !real(dp) JA20,JA21,JA22,JA23,JA24,JA25,JA26,JA27
 !real(dp) JA30,JA31,JA32,JA33,JA34,JA35,JA36
 !real(dp) JA40,JA41,JA42,JA43,JA44,JA45,JA50,JA51,JA52,JA53,JA54
 !real(dp) JA60,JA61,JA62,JA63,JA70,JA71,JA72,JA80,JA81,JA90
 
 real(dp),parameter :: J100=1.e0_dp/3.e0_dp
 
 real(dp),parameter :: J200=1.e0_dp/10.e0_dp
 real(dp),parameter :: J201=2.e0_dp/10.e0_dp
 real(dp),parameter :: J210=1.e0_dp/10.e0_dp
 
 real(dp),parameter :: J300=3.e0_dp/56.e0_dp
 real(dp),parameter :: J301=4.e0_dp/56.e0_dp
 real(dp),parameter :: J302=8.e0_dp/56.e0_dp
 real(dp),parameter :: J310=2.e0_dp/56.e0_dp
 real(dp),parameter :: J311=4.e0_dp/56.e0_dp
 real(dp),parameter :: J320=3.e0_dp/56.e0_dp
 
 real(dp),parameter :: J400=5.e0_dp/144.e0_dp
 real(dp),parameter :: J401=6.e0_dp/144.e0_dp
 real(dp),parameter :: J402=8.e0_dp/144.e0_dp
 real(dp),parameter :: J403=16.e0_dp/144.e0_dp
 real(dp),parameter :: J410=3.e0_dp/144.e0_dp
 real(dp),parameter :: J411=4.e0_dp/144.e0_dp
 real(dp),parameter :: J412=8.e0_dp/144.e0_dp
 real(dp),parameter :: J420=3.e0_dp/144.e0_dp
 real(dp),parameter :: J421=6.e0_dp/144.e0_dp
 real(dp),parameter :: J430=5.e0_dp/144.e0_dp
 
 real(dp),parameter :: J500=35.e0_dp/1408.e0_dp
 real(dp),parameter :: J501=40.e0_dp/1408.e0_dp
 real(dp),parameter :: J502=48.e0_dp/1408.e0_dp
 real(dp),parameter :: J503=64.e0_dp/1408.e0_dp
 real(dp),parameter :: J504=128.e0_dp/1408.e0_dp
 real(dp),parameter :: J510=20.e0_dp/1408.e0_dp
 real(dp),parameter :: J511=24.e0_dp/1408.e0_dp
 real(dp),parameter :: J512=32.e0_dp/1408.e0_dp
 real(dp),parameter :: J513=64.e0_dp/1408.e0_dp
 real(dp),parameter :: J520=18.e0_dp/1408.e0_dp
 real(dp),parameter :: J521=24.e0_dp/1408.e0_dp
 real(dp),parameter :: J522=48.e0_dp/1408.e0_dp
 real(dp),parameter :: J530=20.e0_dp/1408.e0_dp
 real(dp),parameter :: J531=40.e0_dp/1408.e0_dp
 real(dp),parameter :: J540=35.e0_dp/1408.e0_dp
 
 real(dp),parameter :: J600=63.e0_dp/3328.e0_dp
 real(dp),parameter :: J601=70.e0_dp/3328.e0_dp
 real(dp),parameter :: J602=80.e0_dp/3328.e0_dp
 real(dp),parameter :: J603=96.e0_dp/3328.e0_dp
 real(dp),parameter :: J604=128.e0_dp/3328.e0_dp
 real(dp),parameter :: J605=256.e0_dp/3328.e0_dp
 real(dp),parameter :: J610=35.e0_dp/3328.e0_dp
 real(dp),parameter :: J611=40.e0_dp/3328.e0_dp
 real(dp),parameter :: J612=48.e0_dp/3328.e0_dp
 real(dp),parameter :: J613=64.e0_dp/3328.e0_dp
 real(dp),parameter :: J614=128.e0_dp/3328.e0_dp
 real(dp),parameter :: J620=30.e0_dp/3328.e0_dp
 real(dp),parameter :: J621=36.e0_dp/3328.e0_dp
 real(dp),parameter :: J622=48.e0_dp/3328.e0_dp
 real(dp),parameter :: J623=96.e0_dp/3328.e0_dp
 real(dp),parameter :: J630=30.e0_dp/3328.e0_dp
 real(dp),parameter :: J631=40.e0_dp/3328.e0_dp
 real(dp),parameter :: J632=80.e0_dp/3328.e0_dp
 real(dp),parameter :: J640=35.e0_dp/3328.e0_dp
 real(dp),parameter :: J641=70.e0_dp/3328.e0_dp
 real(dp),parameter :: J650=63.e0_dp/3328.e0_dp
 
 real(dp),parameter :: J700=231.e0_dp/15360.e0_dp
 real(dp),parameter :: J701=252.e0_dp/15360.e0_dp
 real(dp),parameter :: J702=280.e0_dp/15360.e0_dp
 real(dp),parameter :: J703=320.e0_dp/15360.e0_dp
 real(dp),parameter :: J704=384.e0_dp/15360.e0_dp
 real(dp),parameter :: J705=512.e0_dp/15360.e0_dp
 real(dp),parameter :: J706=1024.e0_dp/15360.e0_dp
 real(dp),parameter :: J710=126.e0_dp/15360.e0_dp
 real(dp),parameter :: J711=140.e0_dp/15360.e0_dp
 real(dp),parameter :: J712=160.e0_dp/15360.e0_dp
 real(dp),parameter :: J713=192.e0_dp/15360.e0_dp
 real(dp),parameter :: J714=256.e0_dp/15360.e0_dp
 real(dp),parameter :: J715=512.e0_dp/15360.e0_dp
 real(dp),parameter :: J720=105.e0_dp/15360.e0_dp
 real(dp),parameter :: J721=120.e0_dp/15360.e0_dp
 real(dp),parameter :: J722=144.e0_dp/15360.e0_dp
 real(dp),parameter :: J723=192.e0_dp/15360.e0_dp
 real(dp),parameter :: J724=384.e0_dp/15360.e0_dp
 real(dp),parameter :: J730=100.e0_dp/15360.e0_dp
 real(dp),parameter :: J731=120.e0_dp/15360.e0_dp
 real(dp),parameter :: J732=160.e0_dp/15360.e0_dp
 real(dp),parameter :: J733=320.e0_dp/15360.e0_dp
 real(dp),parameter :: J740=105.e0_dp/15360.e0_dp
 real(dp),parameter :: J741=140.e0_dp/15360.e0_dp
 real(dp),parameter :: J742=280.e0_dp/15360.e0_dp
 real(dp),parameter :: J750=126.e0_dp/15360.e0_dp
 real(dp),parameter :: J751=252.e0_dp/15360.e0_dp
 real(dp),parameter :: J760=231.e0_dp/15360.e0_dp
 
 real(dp),parameter :: J800=429.e0_dp/34816.e0_dp
 real(dp),parameter :: J801=462.e0_dp/34816.e0_dp
 real(dp),parameter :: J802=504.e0_dp/34816.e0_dp
 real(dp),parameter :: J803=560.e0_dp/34816.e0_dp
 real(dp),parameter :: J804=640.e0_dp/34816.e0_dp
 real(dp),parameter :: J805=768.e0_dp/34816.e0_dp
 real(dp),parameter :: J806=1024.e0_dp/34816.e0_dp
 real(dp),parameter :: J807=2048.e0_dp/34816.e0_dp
 real(dp),parameter :: J810=231.e0_dp/34816.e0_dp
 real(dp),parameter :: J811=252.e0_dp/34816.e0_dp
 real(dp),parameter :: J812=280.e0_dp/34816.e0_dp
 real(dp),parameter :: J813=320.e0_dp/34816.e0_dp
 real(dp),parameter :: J814=384.e0_dp/34816.e0_dp
 real(dp),parameter :: J815=512.e0_dp/34816.e0_dp
 real(dp),parameter :: J816=1024.e0_dp/34816.e0_dp
 real(dp),parameter :: J820=189.e0_dp/34816.e0_dp
 real(dp),parameter :: J821=210.e0_dp/34816.e0_dp
 real(dp),parameter :: J822=240.e0_dp/34816.e0_dp
 real(dp),parameter :: J823=288.e0_dp/34816.e0_dp
 real(dp),parameter :: J824=284.e0_dp/34816.e0_dp
 real(dp),parameter :: J825=768.e0_dp/34816.e0_dp
 real(dp),parameter :: J830=175.e0_dp/34816.e0_dp
 real(dp),parameter :: J831=200.e0_dp/34816.e0_dp
 real(dp),parameter :: J832=240.e0_dp/34816.e0_dp
 real(dp),parameter :: J833=320.e0_dp/34816.e0_dp
 real(dp),parameter :: J834=640.e0_dp/34816.e0_dp
 real(dp),parameter :: J840=175.e0_dp/34816.e0_dp
 real(dp),parameter :: J841=210.e0_dp/34816.e0_dp
 real(dp),parameter :: J842=280.e0_dp/34816.e0_dp
 real(dp),parameter :: J843=560.e0_dp/34816.e0_dp
 real(dp),parameter :: J850=189.e0_dp/34816.e0_dp
 real(dp),parameter :: J851=252.e0_dp/34816.e0_dp
 real(dp),parameter :: J852=504.e0_dp/34816.e0_dp
 real(dp),parameter :: J860=231.e0_dp/34816.e0_dp
 real(dp),parameter :: J861=462.e0_dp/34816.e0_dp
 real(dp),parameter :: J870=429.e0_dp/34816.e0_dp
 
 real(dp),parameter :: J900=6435.e0_dp/622592.e0_dp
 real(dp),parameter :: J901=6864.e0_dp/622592.e0_dp
 real(dp),parameter :: J902=7392.e0_dp/622592.e0_dp
 real(dp),parameter :: J903=8064.e0_dp/622592.e0_dp
 real(dp),parameter :: J904=8960.e0_dp/622592.e0_dp
 real(dp),parameter :: J905=10240.e0_dp/622592.e0_dp
 real(dp),parameter :: J906=12288.e0_dp/622592.e0_dp
 real(dp),parameter :: J907=16384.e0_dp/622592.e0_dp
 real(dp),parameter :: J908=32768.e0_dp/622592.e0_dp
 real(dp),parameter :: J910=3432.e0_dp/622592.e0_dp
 real(dp),parameter :: J911=3696.e0_dp/622592.e0_dp
 real(dp),parameter :: J912=4032.e0_dp/622592.e0_dp
 real(dp),parameter :: J913=4480.e0_dp/622592.e0_dp
 real(dp),parameter :: J914=5120.e0_dp/622592.e0_dp
 real(dp),parameter :: J915=6144.e0_dp/622592.e0_dp
 real(dp),parameter :: J916=8192.e0_dp/622592.e0_dp
 real(dp),parameter :: J917=16384.e0_dp/622592.e0_dp
 real(dp),parameter :: J920=2772.e0_dp/622592.e0_dp
 real(dp),parameter :: J921=3024.e0_dp/622592.e0_dp
 real(dp),parameter :: J922=3360.e0_dp/622592.e0_dp
 real(dp),parameter :: J923=3840.e0_dp/622592.e0_dp
 real(dp),parameter :: J924=4608.e0_dp/622592.e0_dp
 real(dp),parameter :: J925=6144.e0_dp/622592.e0_dp
 real(dp),parameter :: J926=12288.e0_dp/622592.e0_dp
 real(dp),parameter :: J930=2520.e0_dp/622592.e0_dp
 real(dp),parameter :: J931=2800.e0_dp/622592.e0_dp
 real(dp),parameter :: J932=3200.e0_dp/622592.e0_dp
 real(dp),parameter :: J933=3840.e0_dp/622592.e0_dp
 real(dp),parameter :: J934=5120.e0_dp/622592.e0_dp
 real(dp),parameter :: J935=10240.e0_dp/622592.e0_dp
 real(dp),parameter :: J940=2450.e0_dp/622592.e0_dp
 real(dp),parameter :: J941=2800.e0_dp/622592.e0_dp
 real(dp),parameter :: J942=3360.e0_dp/622592.e0_dp
 real(dp),parameter :: J943=4480.e0_dp/622592.e0_dp
 real(dp),parameter :: J944=8960.e0_dp/622592.e0_dp
 real(dp),parameter :: J950=2520.e0_dp/622592.e0_dp
 real(dp),parameter :: J951=3024.e0_dp/622592.e0_dp
 real(dp),parameter :: J952=4032.e0_dp/622592.e0_dp
 real(dp),parameter :: J953=8064.e0_dp/622592.e0_dp
 real(dp),parameter :: J960=2772.e0_dp/622592.e0_dp
 real(dp),parameter :: J961=3696.e0_dp/622592.e0_dp
 real(dp),parameter :: J962=7392.e0_dp/622592.e0_dp
 real(dp),parameter :: J970=3432.e0_dp/622592.e0_dp
 real(dp),parameter :: J971=6864.e0_dp/622592.e0_dp
 real(dp),parameter :: J980=6435.e0_dp/622592.e0_dp
 
 real(dp),parameter :: JA00=12155.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA01=12870.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA02=13728.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA03=14784.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA04=16128.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA05=17920.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA06=20480.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA07=24576.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA08=32768.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA09=65536.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA10=6435.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA11=6864.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA12=7392.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA13=8064.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA14=8960.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA15=10240.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA16=12288.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA17=16384.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA18=32768.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA20=5148.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA21=5544.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA22=6048.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA23=6720.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA24=7680.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA25=9216.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA26=12288.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA27=24576.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA30=4620.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA31=5040.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA32=5600.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA33=6400.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA34=7680.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA35=10240.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA36=20480.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA40=4410.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA41=4900.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA42=5600.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA43=6720.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA44=8960.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA45=17920.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA50=4410.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA51=5040.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA52=6048.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA53=8064.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA54=16128.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA60=4620.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA61=5544.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA62=7392.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA63=14784.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA70=5148.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA71=6864.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA72=13728.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA80=6435.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA81=12870.e0_dp/1376256.e0_dp
 real(dp),parameter :: JA90=12155.e0_dp/1376256.e0_dp
 
 ! write(*,"(a20,1p3e10.2)") "(serj) y,n,m=",y,n,m
 
 J1=J100
 J2=J200+n*J201+m*J210
 !if(y.lt.9.2925066e-09_dp) then
 !    serj=y*(J1+y*J2)
 ! write(*,"(a20,1pe10.2)") "(serj) J2",serj
 !    return
 !endif
 
 J3=J300+n*(J301+n*J302)+m*(J310+n*J311+m*J320)
 !if(y.lt.4.3667383e-06_dp) then
 !    serj=y*(J1+y*(J2+y*J3))
 ! write(*,"(a20,1pe10.2)") "(serj) J3",serj
 !    return
 !endif
 
 J4=J400+n*(J401+n*(J402+n*J403)) &
     +m*(J410+n*(J411+n*J412)+m*(J420+n*J421+m*J430))
 !if(y.le.9.4990006e-05_dp) then
 !    serj=y*(J1+y*(J2+y*(J3+y*J4)))
 ! write(*,"(a20,1pe10.2)") "(serj) J4",serj
 !    return
 !endif
 
 J5=J500+n*(J501+n*(J502+n*(J503+n*J504))) &
     +m*(J510+n*(J511+n*(J512+n*J513)) &
     +m*(J520+n*(J521+n*J522)+m*(J530+n*J531+m*J540)))
 if(y.le.6.0369310e-04_dp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*J5))))
 ! write(*,"(a20,1pe10.2)") "(serj) J5",serj
     return
 endif
 
 J6=J600+n*(J601+n*(J602+n*(J603+n*(J604+n*J605)))) &
     +m*(J610+n*(J611+n*(J612+n*(J613+n*J614))) &
     +m*(J620+n*(J621+n*(J622+n*J623)) &
     +m*(J630+n*(J631+n*J632)+m*(J640+n*J641+m*J650))))
 if(y.le.2.0727505e-03_dp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*J6)))))
 ! write(*,"(a20,1pe10.2)") "(serj) J6",serj
     return
 endif
 
 J7=J700+n*(J701+n*(J702+n*(J703+n*(J704+n*(J705+n*J706))))) &
     +m*(J710+n*(J711+n*(J712+n*(J713+n*(J714+n*J715)))) &
     +m*(J720+n*(J721+n*(J722+n*(J723+n*J724))) &
     +m*(J730+n*(J731+n*(J732+n*J733)) &
     +m*(J740+n*(J741+n*J742)+m*(J750+n*J751+m*J760)))))
 if(y.le.5.0047026e-03_dp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*J7))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J7",serj
    return
 endif
 
 J8=J800+n*(J801+n*(J802+n*(J803+n*(J804+n*(J805+n*(J806+n*J807)))))) &
     +m*(J810+n*(J811+n*(J812+n*(J813+n*(J814+n*(J815+n*J816))))) &
     +m*(J820+n*(J821+n*(J822+n*(J823+n*(J824+n*J825)))) &
     +m*(J830+n*(J831+n*(J832+n*(J833+n*J834))) &
     +m*(J840+n*(J841+n*(J842+n*J843)) &
     +m*(J850+n*(J851+n*J852)+m*(J860+n*J861+m*J870))))))
 if(y.le.9.6961652e-03_dp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*J8)))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J8",serj
     return
 endif
 
 J9=J900+n*(J901+n*(J902+n*(J903+n*(J904+n*(J905+n*(J906+n*(J907+n*J908))))))) &
     +m*(J910+n*(J911+n*(J912+n*(J913+n*(J914+n*(J915+n*(J916+n*J917)))))) &
     +m*(J920+n*(J921+n*(J922+n*(J923+n*(J924+n*(J925+n*J926))))) &
     +m*(J930+n*(J931+n*(J932+n*(J933+n*(J934+n*J935)))) &
     +m*(J940+n*(J941+n*(J942+n*(J943+n*J944))) &
     +m*(J950+n*(J951+n*(J952+n*J953)) &
     +m*(J960+n*(J961+n*J962)+m*(J970+n*J971+m*J980)))))))
 if(y.le.1.6220210e-02_dp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*(J8+y*J9))))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J9",serj
     return
 endif
 
 JA=JA00+n*(JA01+n*(JA02+n*(JA03+n*(JA04+n*(JA05+n*(JA06+n*(JA07+n*(JA08+n*JA09)))))))) &
     +m*(JA10+n*(JA11+n*(JA12+n*(JA13+n*(JA14+n*(JA15+n*(JA16+n*(JA17+n*JA18))))))) &
     +m*(JA20+n*(JA21+n*(JA22+n*(JA23+n*(JA24+n*(JA25+n*(JA26+n*JA27)))))) &
     +m*(JA30+n*(JA31+n*(JA32+n*(JA33+n*(JA34+n*(JA35+n*JA36))))) &
     +m*(JA40+n*(JA41+n*(JA42+n*(JA43+n*(JA44+n*JA45)))) &
     +m*(JA50+n*(JA51+n*(JA52+n*(JA53+n*JA54))) &
     +m*(JA60+n*(JA61+n*(JA62+n*JA63)) &
     +m*(JA70+n*(JA71+n*JA72)+m*(JA80+n*JA81+m*JA90))))))))
 !if(y.lt.2.4482875d-02) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*(J8+y*(J9+y*JA)))))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J10",serj
     return
 !endif
 
 end
 !---------------------------------------------------------------------------
 elemental real(dp) function uatan(t,h)
 !
 ! Universal arctangent function
 !
 ! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !   Precise and fast computation of a general incomplete elliptic integral
 !   of third kind by half and double argument transformations
 !
 real(dp),intent(in) :: t,h !!SJT
 real(dp) :: z,y,x !!SJT
 real(dp) a,r,ri,hold,rold,riold
 !real(dp) A3,A5,A7,A9,A11,A13,A15,A17,A19,A21,A23,A25
 !data hold/1.d0/, rold/1.d0/,riold/1.d0/ !!SJT: disable for thread safety
 !save hold,rold,riold !!SJT: disable for thread safety
 real(dp),parameter :: A3=1.e0_dp/3.e0_dp
 real(dp),parameter :: A5=1.e0_dp/5.e0_dp
 real(dp),parameter :: A7=1.e0_dp/7.e0_dp
 real(dp),parameter :: A9=1.e0_dp/9.e0_dp
 real(dp),parameter :: A11=1.e0_dp/11.e0_dp
 real(dp),parameter :: A13=1.e0_dp/13.e0_dp
 real(dp),parameter :: A15=1.e0_dp/15.e0_dp
 real(dp),parameter :: A17=1.e0_dp/17.e0_dp
 real(dp),parameter :: A19=1.e0_dp/19.e0_dp
 real(dp),parameter :: A21=1.e0_dp/21.e0_dp
 real(dp),parameter :: A23=1.e0_dp/23.e0_dp
 real(dp),parameter :: A25=1.e0_dp/25.e0_dp
 !!!!!SJT: replace data statment with manual value assignment for thread safety
 hold=1.e0_dp; rold=1.e0_dp; riold=1.e0_dp
 !!!!!SJT: end replace data statment with manual value assignment for thread safety
 !
 ! write(*,*) "(uatan) t,h=",t,h
 !
 z=-h*t*t
 a=abs(z)
 !
 ! write(*,*) "(uatan) z=",z
 if(a.lt.3.3306691e-16_dp) then
     uatan=t
 ! write(*,*) "0"
 elseif(a.lt.2.3560805e-08_dp) then
     uatan=t*(1.e0_dp+z*A3)
 ! write(*,*) "1"
 elseif(a.lt.9.1939631e-06_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*A5))
 ! write(*,*) "2"
 elseif(a.lt.1.7779240e-04_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*A7)))
 ! write(*,*) "3"
 elseif(a.lt.1.0407839e-03_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*(A7+z*A9))))
 ! write(*,*) "4"
 elseif(a.lt.3.3616998e-03_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*(A7+z*(A9+z*A11)))))
 ! write(*,*) "5"
 elseif(a.lt.7.7408014e-03_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*A13))))))
 ! write(*,*) "6"
 elseif(a.lt.1.4437181e-02_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*A15)))))))
 ! write(*,*) "7"
 elseif(a.lt.2.3407312e-02_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*A17))))))))
 ! write(*,*) "8"
 elseif(a.lt.3.4416203e-02_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*A19)))))))))
 ! write(*,*) "9"
 elseif(z.lt.0.e0_dp) then
     if(abs(h-hold).lt.1.e-16_dp) then
         r=rold; ri=riold
     else
         r=sqrt(h); ri=1.e0_dp/r; hold=h; rold=r; riold=ri
     endif
     uatan=atan(r*t)*ri
 ! write(*,*) "atan"
 elseif(a.lt.4.7138547e-02_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*A21))))))))))
 ! write(*,*) "A"
 elseif(a.lt.6.1227405e-02_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*(A21+z*A23)))))))))))
 ! write(*,*) "B"
 elseif(a.lt.7.6353468e-02_dp) then
     uatan=t*(1.e0_dp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*(A21+z*(A23+z*A25))))))))))))
 ! write(*,*) "C"
 else
     if(abs(h-hold).lt.1.e-16_dp) then
         r=rold; ri=riold
     else
         r=sqrt(-h); ri=1.e0_dp/r; hold=h; rold=r; riold=ri
     endif
     y=r*t
     x=log((1.e0_dp+y)/(1.e0_dp-y))*0.5e0_dp
     uatan=x*ri
 ! write(*,*) "hyper"
 endif
 !
 ! write(*,*) "(uatan) uatan=",uatan
 !
 return
 end
 

! *************************************** SJT: Extra Precision Routines Start ***************************************

 elemental subroutine elbdj2_qp(phi,phic,n,mc_qp,b,d,j) !!SJT
 !subroutine elbdj2(phi,phic,n,mc,b,d,j) !!SJT: OG
 !
 !     Simultaneous computation of associate elliptic integrals
 !     of third kind, B(phi|m), D(phi|m), and J(phi,n|m)
 !
 !     Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !       Precise and fast computation of a general incomplete elliptic integral
 !       of third kind by half and double argument transformations
 !
 !     Modified Version by employing "elsbdj" only
 !
 !     Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
 !
 !     Used subprograms: celbd,celbdj,elsbdj,serbd,serj,uatan
 !
 !     Inputs: phi  = argument                0 <= phi  <= PI/2
 !             phic = complementar argument   0 <= phic <= PI/2
 !             n    = characteristic          0 <= n    <= 1
 !             mc   = complementary parameter 0 <= mc   <= 1
 !
 !     Outputs: b, d, j
 !
 !     CAUTION: phi and phic must satisfy condition, phi + phic = PI/2
 !
 real(qp),intent(in) :: phi,phic,n !!SJT
 real(qp),intent(out) :: b,d,j !!SJT
 real(qp) mc,m,nc,h,c,x,d2,z,sz,t,v,t2 !!SJT
 real(qp) bc,dc,jc
 real(qp),intent(in) :: mc_qp !!SJT
 !
 if(phi.lt.1.249e0_qp) then   ! Modified: Pascal Leroy's suggestion 2019/08/12
     call elsbdj_qp(sin(phi),n,mc_qp,b,d,j) !!SJT
 else
     !m=1.d0-mc !!SJT: original
     m=1.0_qp-mc_qp; mc=mc_qp !!SJT
     nc=1.e0_qp-n
     h=n*nc*(n-m)
     c=sin(phic)
     x=c*c
     d2=mc+m*x
     if(x.lt.0.9e0_qp*d2) then
         z=c/sqrt(d2)
         call elsbdj_qp(z,n,mc_qp,b,d,j) !!SJT
         call celbdj_qp(nc,mc_qp,bc,dc,jc) !!SJT
                 sz=z*sqrt(1.e0_qp-x)
                 t=sz/nc
                 b=bc-(b-sz)
                 d=dc-(d+sz)
                 j=jc-(j+uatan_qp(t,h))
     else
                 v=mc*(1.e0_qp-x)
         if(v.lt.x*d2) then
             call elsbdj_qp(sqrt(1.e0_qp-c*c),n,mc_qp,b,d,j) !!SJT
         else
             t2=(1.e0_qp-x)/d2
             call elsbdj_qp(sqrt(1.e0_qp-mc*t2),n,mc_qp,b,d,j) !!SJT
             call celbdj_qp(nc,mc_qp,bc,dc,jc) !!SJT
                     sz=c*sqrt(t2)
                     t=sz/nc
                     b=bc-(b-sz)
                     d=dc-(d+sz)
                     j=jc-(j+uatan_qp(t,h))
         endif
     endif
 endif
 return
 end

 elemental subroutine elsbdj_qp(s0,n,mc_qp,b,d,j) !!SJT
 !subroutine elsbdj(s0,n,mc,b,d,j) !!SJT: OG
 !
 ! Simultaneous computation of associate elliptic integrals
 ! of third kind, B(phi|m), D(phi|m), and J(phi,n|m)
 ! by using the half/double argument transformation of sn functions
 !
 ! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !   Precise and fast computation of a general incomplete elliptic integral
 !   of third kind by half and double argument transformations
 !
 real(qp),intent(in) :: s0,n !!SJT
 real(qp),intent(out) :: b,d,j !!SJT
 real(qp) m,h,del,s,y,c,sy,t !!SJT
 real(qp),intent(in) :: mc_qp !!SJT
 integer(isp),parameter :: NMAX=40_isp ! SJT: added variable max size (original value is 10)
 !real(qp) yy(11),ss(11),cd(11) ! OG
 real(qp) yy(NMAX+1_isp),ss(NMAX+1_isp),cd(NMAX+1_isp)
 integer(isp) i,k
 
 ! write(*,*) "(elsj) s0,n,mc=",s0,n,mc
 
 !m=1.d0-mc !!SJT: original
 m=1.0_qp-mc_qp 
 h=n*(1._qp-n)*(n-m)
 ! SJT: note that the J9 term was used originally (enough for double precision accuracy)
 ! SJT: smaller del values allow more terms to be retained, improving precision
 !del=0.04094_qp-0.00652_qp*m
 !del=0.03429_qp
 !del=0.02448_qp  ! JA
 !del=0.01622_qp   ! J9
 !del=0.00969_qp  ! J8
 !del=0.00500_qp  ! J7
 !del=2.073e-3_qp    ! J6
 del=6.037e-4_qp    ! J5
 
 s=s0
 y=s*s
 if(y.lt.del) then
         call serbd_qp(y,m,b,d)
         b=s*b
         d=s*y*d
         j=s*serj_qp(y,n,m)
 ! write(*,"(a20,1p3e10.2)") "(elsbdj) b,d,j=",b,d,j
         return
 endif
 yy(1_isp)=y
 ss(1_isp)=s
 ! write(*,"(a20,i10,1pe10.2)") "(elsbdj) i,y=",1,y
 !do i=1,10 ! OG
 do i=1_isp,NMAX
     c=sqrt(1._qp-y)
     d=sqrt(1._qp-m*y)
     y=y/((1._qp+c)*(1._qp+d))
         yy(i+1_isp)=y
         ss(i+1_isp)=sqrt(y)
         cd(i)=c*d
 ! write(*,"(a20,i10,1pe10.2)") "(elsbdj) i,y=",i+1,y
         if(y.lt.del) then
                 goto 1
         endif
 enddo
 !write(*,*) "(elsbdj) too many iterations: s0,m=",s0,m SJT: does not occur in practice
 1 continue
 ! write(*,"(a20,i10,1pe10.2)") "(elsbdj) i,y=",i+1,y
 call serbd_qp(y,m,b,d)
 b=ss(i+1_isp)*b
 d=ss(i+1_isp)*y*d
 j=ss(i+1_isp)*serj_qp(y,n,m)
 do k=i,1_isp,-1_isp
         sy=ss(k)*yy(k+1_isp)
         t=sy/(1._qp-n*(yy(k)-yy(k+1_isp)*cd(k)))
         b=2._qp*b-sy
         d=d+(d+sy)
         j=j+(j+uatan_qp(t,h))
 enddo
 ! write(*,"(a20,1p3e10.2)") "(elsbdj) b,d,j=",b,d,j
 return
 end

 elemental subroutine celbdj_qp(nc0,mc0_qp,celb,celd,celj)
 !subroutine celbdj(nc0,mc0,celb,celd,celj) !!SJT: OG
 !
 ! Simultaneous computation of associate complete elliptic integrals
 ! of third kind, B(m), D(m), and J(n|m)
 !
 ! Reference: Fukushima, T, (2013) J. Comp. Appl. Math., 253, 142-157
 !     Fast computation of a general complete elliptic integral
 !     of third kind by half and double argument transformations
 !
 ! Restricted Version with Assumptions: 0 < mc, 0 < nc
 !
 ! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
 !
 real(qp),intent(in) :: nc0
 real(qp),intent(out) :: celb,celd,celj !!SJT
 logical flag
 !integer(isp) IMAX,i,is,ie !!SJT: OG
 integer(isp) i,is,ie
 !parameter (IMAX=40) ! OG
 integer(isp),parameter :: IMAX=40_isp  ! SJT: adjustement in max # of iterations for quadruple precision
 real(qp) y(0_isp:IMAX),x(0_isp:IMAX),c(0_isp:IMAX),d(0_isp:IMAX),a(0_isp:IMAX)
 real(qp) mc,mc0,nc,m,n,celk,yi,ye,dj,m1,kc0,temp !!SJT
 real(qp),intent(in) :: mc0_qp !!SJT
 real(qp) :: mc_qp !!SJT
 !real(qp) PIHALF,EPS,THIRD
 real(qp),parameter :: PIHALF=1.570796326794896619231321691639751442099_qp
 !parameter (EPS=1.11e-16_qp) ! SJT: OG value
 real(qp),parameter :: EPS=1.e-34_qp ! SJT: adjusted for quadruple precision (roughly half of machine epsilon)
 real(qp),parameter :: THIRD=1.e0_qp/3.e0_qp
 real(qp) B(IMAX)
 !logical first /.TRUE./  !!!SJT: note that logical variable first is .true. at the start of each routine call
 !save B,first !!SJT: disable for thread safety
 logical first !!SJT: initialization cannot be in declaration for thread-safe calculations
 first=.true. !!SJT: initialization outside of declaration
 mc0=mc0_qp !!SJT
 if(mc0.lt.1._qp) then
     !mc=mc0;nc=nc0
     mc_qp=mc0_qp; nc=nc0 !!SJT
 elseif(mc0.gt.1._qp) then
     !mc=1.d0/mc0;nc=nc0*mc !!SJT: OG
     mc_qp=1.0_qp/mc0_qp; nc=nc0*mc_qp !!SJT
 else
     celb=PIHALF*0.5_qp
     celd=PIHALF*0.5_qp
     celj=PIHALF/(nc0+sqrt(nc0))
     return
 endif
 if(first) then
     first=.FALSE.
     do i=1,IMAX
         B(i)=1._qp/real(2_isp*i+1_isp,qp)
     enddo
 endif
 mc=mc_qp !!SJT
 !call celbd(mc,celb,celd) !!SJT: OG
 call celbd_qp(mc_qp,celb,celd) !!SJT
 if(nc.eq.mc) then
     celj=celb/mc
     goto 4  
 endif
 !m=1.d0-mc; n=1.d0-nc !!SJT: original
 m=1.0_qp-mc_qp; n=1._qp-nc !!SJT
 flag=nc.lt.mc.or.(n*nc).gt.(nc-mc)
 if(flag) then
     y(0_isp)=(nc-mc)/(nc*m)
 else
     y(0_isp)=n/m
 endif
 is=0_isp
 if(y(0_isp).gt.0.5_qp) then
     x(0_isp)=1._qp-y(0_isp)
     do i=0_isp,IMAX-1_isp !!SJT
     !do i=0,IMAX !!SJT: OG
         c(i)=sqrt(x(i))
         d(i)=sqrt(mc+m*x(i))
         x(i+1)=(c(i)+d(i))/(1._qp+d(i))
         if(x(i+1).gt.0.5_qp) then
             y(i+1)=y(i)/((1._qp+c(i))*(1._qp+d(i)))
             is=i+1_isp
             goto 1
         endif
         y(i+1_isp)=1._qp-x(i+1_isp)
     enddo
     return
 endif
 1 continue
 do i=is,IMAX-1_isp !!SJT
 !do i=is,IMAX !!SJT: OG
     c(i)=sqrt(1._qp-y(i))
     d(i)=sqrt(1._qp-m*y(i))
     y(i+1)=y(i)/((1._qp+c(i))*(1._qp+d(i)))
     if(abs(y(i+1)).lt.0.325_qp) goto 2
 enddo
 return
 2 continue
 ie=i+1_isp
 ye=y(ie)
 celk=celb+celd
 a(0_isp)=celd
 celj=a(0_isp)
 yi=ye
 a(1_isp)=((1._qp+2._qp*m)*celd-celb)*THIRD
 dj=a(1_isp)*yi
 i=1_isp
 if(abs(dj).lt.EPS*abs(celj)) goto 3
 celj=celj+dj
 m1=1._qp+m
 do i=2_isp,IMAX
     yi=yi*ye
     a(i)=(1._qp-B(i))*m1*a(i-1_isp)-(1._qp-2._qp*B(i))*m*a(i-2_isp)
     dj=a(i)*yi
     if(abs(dj).lt.EPS*abs(celj)) goto 3
     celj=celj+dj
 enddo
 return
 3 continue
 do i=ie-1_isp,0_isp,-1_isp
     celj=(2._qp*(c(i)+d(i))*celj-y(i)*celk)/(c(i)*d(i)*(1._qp+c(i))*(1._qp+d(i)))
 enddo
 if(flag) then
     celj=(nc*celk-mc*celj)/(nc*nc)
 endif
 4 continue
 if(mc0.gt.1._qp) then
     kc0=sqrt(mc0)
     temp=celb
     celb=celd/kc0
     celd=temp/kc0
     celj=celj/(mc0*kc0)
 endif
 return
 end

 elemental subroutine celbd_qp(mc_qp,elb,eld) !!SJT
 !subroutine celbd(mc,elb,eld) !!SJT: OG
 !
 ! Simultaneous computation of associate complete elliptic integrals
 ! of second kind, B(m) and D(m)
 !
 ! Reference: Fukushima, T (2011) Math. Comp., 80, 1725-1743
 !   Precise and fast computation of general complete elliptic integral
 !   of second kind
 !
 ! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
 !
 real(qp) mc,elk !!SJT
 real(qp),intent(out) :: elb,eld !!SJT
 real(qp) m,mx,kkc,nome
 real(qp),intent(in) :: mc_qp !!SJT
 !
 real(qp),parameter :: EPS=1.e-34_qp  !1.11e-16_qp ! SJT: adjusted for quadruple precision (~half of machine epsilon) 
 !real(qp) PIQ,PIHALF!,PI,PIINV !!SJT: removed unused parameters
 real(qp),parameter :: PIQ=0.78539816339744830961566084581988_qp
 real(qp),parameter :: PIHALF=1.5707963267948966192313216916398_qp
 !parameter (PI=3.1415926535897932384626433832795d0)
 !parameter (PIINV=0.31830988618379067153776752674503d0)
 real(qp) mcold,elbold,eldold
 !save mcold,elbold,eldold !!SJT: disable for thread safety -- variables are initialized at start of each celbd call
 !real(qp) Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16
 real(qp),parameter :: Q1=1._qp/16._qp,Q2=1._qp/32._qp,Q3=21._qp/1024._qp
 real(qp),parameter :: Q4=31._qp/2048._qp,Q5=6257._qp/524288._qp
 real(qp),parameter :: Q6=10293._qp/1048576._qp,Q7=279025._qp/33554432._qp
 real(qp),parameter :: Q8=483127._qp/67108864._qp
 real(qp),parameter :: Q9=435506703._qp/68719476736._qp
 real(qp),parameter :: Q10=776957575._qp/137438953472._qp
 real(qp),parameter :: Q11=22417045555._qp/4398046511104._qp
 real(qp),parameter :: Q12=40784671953._qp/8796093022208._qp
 real(qp),parameter :: Q13=9569130097211._qp/2251799813685248._qp
 real(qp),parameter :: Q14=17652604545791._qp/4503599627370496._qp
 real(qp),parameter :: Q15=523910972020563._qp/144115188075855872._qp
 real(qp),parameter :: Q16=976501268709949._qp/288230376151711744._qp
 !real(qp) K1,K2,K3,K4,K5,K6,K7
 real(qp),parameter :: K1=1._qp/4._qp
 real(qp),parameter :: K2=9._qp/64._qp
 real(qp),parameter :: K3=25._qp/256._qp
 real(qp),parameter :: K4=1225._qp/16384._qp
 real(qp),parameter :: K5=3969._qp/65536._qp
 real(qp),parameter :: K6=53361._qp/1048576._qp
 real(qp),parameter :: K7=184041._qp/4194304._qp
 !real(qp) B1,B2,B3,B4,B5,B6,B7,B8
 real(qp),parameter :: B1=1._qp/2._qp
 real(qp),parameter :: B2=1._qp/16._qp
 real(qp),parameter :: B3=3._qp/128._qp
 real(qp),parameter :: B4=25._qp/2048._qp
 real(qp),parameter :: B5=245._qp/32768._qp
 real(qp),parameter :: B6=1323._qp/262144._qp
 real(qp),parameter :: B7=7623._qp/2097152._qp
 real(qp),parameter :: B8=184041._qp/67108864._qp
 !real(qp) D1,D2,D3,D4,D5,D6,D7,D8
 real(qp),parameter :: D1=1._qp/2._qp
 real(qp),parameter :: D2=3._qp/16._qp
 real(qp),parameter :: D3=15._qp/128._qp
 real(qp),parameter :: D4=175._qp/2048._qp
 real(qp),parameter :: D5=2205._qp/32768._qp
 real(qp),parameter :: D6=14553._qp/262144._qp
 real(qp),parameter :: D7=99099._qp/2097152._qp
 real(qp),parameter :: D8=2760615._qp/67108864._qp
 real(qp) logq2,dkkc,dddc,dele,delb,elk1
 !logical first/.TRUE./    !!!SJT: note that logical variable first is always .true. at the start of each routine call
 logical first !!SJT: initialization within declaration not allowed due to thread safety
 first=.true.  !!SJT: initialize outside of declaration for thread safety
 if(first) then
     first=.FALSE.
         mcold=1._qp
         elbold=PIQ
         eldold=PIQ
 endif
 !m=1.d0-mc !!SJT: original
 m=1.0_qp-mc_qp; mc=mc_qp !!SJT
 !if(abs(mc-mcold).lt.1.11e-16_qp*mc) then
 if(abs(mc-mcold).lt.EPS*mc) then
         elb=elbold
         eld=eldold
 !elseif(m.lt.1.11e-16_qp) then
 elseif(m.lt.EPS) then
     elb=PIQ
         eld=PIQ
 !elseif(mc.lt.1.11e-16_qp) then
 elseif(mc.lt.EPS) then
     elb=1._qp
         eld=0.3862943611198906188344642429164e0_qp-0.5e0_qp*log(mc)
 elseif(mc.lt.0.1e0_qp) then
     nome=mc*(Q1+mc*(Q2+mc*(Q3+mc*(Q4+mc*(Q5+mc*(Q6 &
         +mc*(Q7+mc*(Q8+mc*(Q9+mc*(Q10+mc*(Q11+mc*(Q12 &
         +mc*(Q13+mc*(Q14+mc*(Q15+mc*Q16))))))))))))))) 
     if(mc.lt.0.01e0_qp) then
         dkkc=mc*(K1+mc*(K2+mc*(K3+mc*(K4+mc*(K5+mc*(K6+mc*K7))))))
         dddc=mc*(D1+mc*(D2+mc*(D3+mc*(D4+mc*(D5+mc*(D6+mc*D7))))))
     else
         mx=mc-0.05e0_qp
 !
 ! (K'-1)/(pi/2)
 !
         dkkc=    0.01286425658832983978282698630501405107893_qp &
             +mx*(0.26483429894479586582278131697637750604652_qp &
             +mx*(0.15647573786069663900214275050014481397750_qp &
             +mx*(0.11426146079748350067910196981167739749361_qp &
             +mx*(0.09202724415743445309239690377424239940545_qp &
             +mx*(0.07843218831801764082998285878311322932444_qp &
             +mx*(0.06935260142642158347117402021639363379689_qp &
             +mx*(0.06293203529021269706312943517695310879457_qp &
             +mx*(0.05821227592779397036582491084172892108196_qp &
             +mx*(0.05464909112091564816652510649708377642504_qp &
             +mx*(0.05191068843704411873477650167894906357568_qp &
             +mx*(0.04978344771840508342564702588639140680363_qp &
             +mx*(0.04812375496807025605361215168677991360500_qp &
             ))))))))))))
 !
 ! (K'-E')/(pi/2)
 !
         dddc=    0.02548395442966088473597712420249483947953_qp &
             +mx*(0.51967384324140471318255255900132590084179_qp &
             +mx*(0.20644951110163173131719312525729037023377_qp &
             +mx*(0.13610952125712137420240739057403788152260_qp &
             +mx*(0.10458014040566978574883406877392984277718_qp &
             +mx*(0.08674612915759188982465635633597382093113_qp &
             +mx*(0.07536380269617058326770965489534014190391_qp &
             +mx*(0.06754544594618781950496091910264174396541_qp &
             +mx*(0.06190939688096410201497509102047998554900_qp &
             +mx*(0.05771071515451786553160533778648705873199_qp &
             +mx*(0.05451217098672207169493767625617704078257_qp &
             +mx*(0.05204028407582600387265992107877094920787_qp &
             +mx*(0.05011532514520838441892567405879742720039_qp &
             ))))))))))))
     endif
         kkc=1._qp+dkkc
         logq2=-0.5_qp*log(nome)
         elk=kkc*logq2
         dele=-dkkc/kkc+logq2*dddc
         elk1=elk-1._qp
         delb=(dele-mc*elk1)/m
         elb=1._qp+delb
         eld=elk1-delb
 elseif(m.le.0.01_qp) then
         elb=PIHALF*(B1+m*(B2+m*(B3+m*(B4+m*(B5+m*(B6+m*(B7+m*B8)))))))
         eld=PIHALF*(D1+m*(D2+m*(D3+m*(D4+m*(D5+m*(D6+m*(D7+m*D8)))))))
 elseif(m.le.0.1_qp) then
         mx=0.95_qp-mc
         elb=     0.790401413584395132310045630540381158921005_qp &
             +mx*(0.102006266220019154892513446364386528537788_qp &
             +mx*(0.039878395558551460860377468871167215878458_qp &
             +mx*(0.021737136375982167333478696987134316809322_qp &
             +mx*(0.013960979767622057852185340153691548520857_qp &
             +mx*(0.009892518822669142478846083436285145400444_qp &
             +mx*(0.007484612400663335676130416571517444936951_qp &
             +mx*(0.005934625664295473695080715589652011420808_qp &
             +mx*(0.004874249053581664096949448689997843978535_qp &
             +mx*(0.004114606930310886136960940893002069423559_qp &
             +mx*(0.003550452989196176932747744728766021440856_qp &
             +mx*(0.003119229959988474753291950759202798352266_qp &
             )))))))))))
         eld=     0.800602040206397047799296975176819811774784_qp &
             +mx*(0.313994477771767756849615832867393028789057_qp &
             +mx*(0.205913118705551954501930953451976374435088_qp &
             +mx*(0.157744346538923994475225014971416837073598_qp &
             +mx*(0.130595077319933091909091103101366509387938_qp &
             +mx*(0.113308474489758568672985167742047066367053_qp &
             +mx*(0.101454199173630195376251916342483192174927_qp &
             +mx*(0.0929187842072974367037702927967784464949434_qp &
             +mx*(0.0865653801481680871714054745336652101162894_qp &
             +mx*(0.0817279846651030135350056216958053404884715_qp &
             +mx*(0.0779906657291070378163237851392095284454654_qp &
             +mx*(0.075080426851268007156477347905308063808848_qp &
             )))))))))))
 elseif(m.le.0.2_qp) then
         mx=0.85_qp-mc
         elb=     0.80102406445284489393880821604009991524037_qp &
             +mx*(0.11069534452963401497502459778015097487115_qp &
             +mx*(0.047348746716993717753569559936346358937777_qp &
             +mx*(0.028484367255041422845322166419447281776162_qp &
             +mx*(0.020277811444003597057721308432225505126013_qp &
             +mx*(0.015965005853099119442287313909177068173564_qp &
             +mx*(0.013441320273553634762716845175446390822633_qp &
             +mx*(0.011871565736951439501853534319081030547931_qp &
             +mx*(0.010868363672485520630005005782151743785248_qp &
             +mx*(0.010231587232710564565903812652581252337697_qp &
             +mx*(0.009849585546666211201566987057592610884309_qp &
             +mx*(0.009656606347153765129943681090056980586986_qp &
             )))))))))))
         eld=     0.834232667811735098431315595374145207701720_qp &
             +mx*(0.360495281619098275577215529302260739976126_qp &
             +mx*(0.262379664114505869328637749459234348287432_qp &
             +mx*(0.223723944518094276386520735054801578584350_qp &
             +mx*(0.206447811775681052682922746753795148394463_qp &
             +mx*(0.199809440876486856438050774316751253389944_qp &
             +mx*(0.199667451603795274869211409350873244844882_qp &
             +mx*(0.204157558868236842039815028663379643303565_qp &
             +mx*(0.212387467960572375038025392458549025660994_qp &
             +mx*(0.223948914061499360356873401571821627069173_qp &
             +mx*(0.238708097425597860161720875806632864507536_qp &
             +mx*(0.256707203545463755643710021815937785120030_qp &
             )))))))))))
 elseif(m.le.0.3_qp) then
         mx=0.75_qp-mc
         elb=     0.81259777291992049322557009456643357559904_qp &
             +mx*(0.12110961794551011284012693733241967660542_qp &
             +mx*(0.057293376831239877456538980381277010644332_qp &
             +mx*(0.038509451602167328057004166642521093142114_qp &
             +mx*(0.030783430301775232744816612250699163538318_qp &
             +mx*(0.027290564934732526869467118496664914274956_qp &
             +mx*(0.025916369289445198731886546557337255438083_qp &
             +mx*(0.025847203343361799141092472018796130324244_qp &
             +mx*(0.026740923539348854616932735567182946385269_qp &
             +mx*(0.028464314554825704963640157657034405579849_qp &
             +mx*(0.030995446237278954096189768338119395563447_qp &
             +mx*(0.034384369179940975864103666824736551261799_qp &
             +mx*(0.038738002072493935952384233588242422046537_qp &
             ))))))))))))
         eld=     0.873152581892675549645633563232643413901757_qp &
             +mx*(0.420622230667770215976919792378536040460605_qp &
             +mx*(0.344231061559450379368201151870166692934830_qp &
             +mx*(0.331133021818721761888662390999045979071436_qp &
             +mx*(0.345277285052808411877017306497954757532251_qp &
             +mx*(0.377945322150393391759797943135325823338761_qp &
             +mx*(0.427378012464553880508348757311170776829930_qp &
             +mx*(0.494671744307822405584118022550673740404732_qp &
             +mx*(0.582685115665646200824237214098764913658889_qp &
             +mx*(0.695799207728083164790111837174250683834359_qp &
             +mx*(0.840018401472533403272555302636558338772258_qp &
             +mx*(1.023268503573606060588689738498395211300483_qp &
             +mx*(1.255859085136282496149035687741403985044122_qp &
             ))))))))))))
 elseif(m.le.0.4_qp) then
         mx=0.65_qp-mc
         elb=     0.8253235579835158949845697805395190063745_qp &
             +mx*(0.1338621160836877898575391383950840569989_qp &
             +mx*(0.0710112935979886745743770664203746758134_qp &
             +mx*(0.0541784774173873762208472152701393154906_qp &
             +mx*(0.0494517449481029932714386586401273353617_qp &
             +mx*(0.0502221962241074764652127892365024315554_qp &
             +mx*(0.0547429131718303528104722303305931350375_qp &
             +mx*(0.0627462579270016992000788492778894700075_qp &
             +mx*(0.0746698810434768864678760362745179321956_qp &
             +mx*(0.0914808451777334717996463421986810092918_qp &
             +mx*(0.1147050921109978235104185800057554574708_qp &
             +mx*(0.1465711325814398757043492181099197917984_qp &
             +mx*(0.1902571373338462844225085057953823854177_qp &
             ))))))))))))
         eld=     0.9190270392420973478848471774160778462738_qp &
             +mx*(0.5010021592882475139767453081737767171354_qp &
             +mx*(0.4688312705664568629356644841691659415972_qp &
             +mx*(0.5177142277764000147059587510833317474467_qp &
             +mx*(0.6208433913173031070711926900889045286988_qp &
             +mx*(0.7823643937868697229213240489900179142670_qp &
             +mx*(1.0191145350761029126165253557593691585239_qp &
             +mx*(1.3593452027484960522212885423056424704073_qp &
             +mx*(1.8457173023588279422916645725184952058635_qp &
             +mx*(2.5410717031539207287662105618152273788399_qp &
             +mx*(3.5374046552080413366422791595672470037341_qp &
             +mx*(4.9692960029774259303491034652093672488707_qp &
             +mx*(7.0338228700300311264031522795337352226926_qp &
             +mx*(10.020043225034471401553194050933390974016_qp &
             )))))))))))))
 elseif(m.le.0.5_qp) then
         mx=0.55_qp-mc
         elb=     0.8394795702706129706783934654948360410325_qp &
             +mx*(0.1499164403063963359478614453083470750543_qp &
             +mx*(0.0908319358194288345999005586556105610069_qp &
             +mx*(0.0803470334833417864262134081954987019902_qp &
             +mx*(0.0856384405004704542717663971835424473169_qp &
             +mx*(0.1019547259329903716766105911448528069506_qp &
             +mx*(0.1305748115336160150072309911623351523284_qp &
             +mx*(0.1761050763588499277679704537732929242811_qp &
             +mx*(0.2468351644029554468698889593583314853486_qp &
             +mx*(0.3564244768677188553323196975301769697977_qp &
             +mx*(0.5270025622301027434418321205779314762241_qp &
             +mx*(0.7943896342593047502260866957039427731776_qp &
             +mx*(1.2167625324297180206378753787253096783993_qp &
             ))))))))))))
         eld=     0.9744043665463696730314687662723484085813_qp &
             +mx*(0.6132468053941609101234053415051402349752_qp &
             +mx*(0.6710966695021669963502789954058993004082_qp &
             +mx*(0.8707276201850861403618528872292437242726_qp &
             +mx*(1.2295422312026907609906452348037196571302_qp &
             +mx*(1.8266059675444205694817638548699906990301_qp &
             +mx*(2.8069345309977627400322167438821024032409_qp &
             +mx*(4.4187893290840281339827573139793805587268_qp &
             +mx*(7.0832360574787653249799018590860687062869_qp &
             +mx*(11.515088120557582942290563338274745712174_qp &
             +mx*(18.931511185999274639516732819605594455165_qp &
             +mx*(31.411996938204963878089048091424028309798_qp &
             +mx*(52.520729454575828537934780076768577185134_qp &
             +mx*(88.384854735065298062125622417251073520996_qp &
             +mx*(149.56637449398047835236703116483092644714_qp &
             +mx*(254.31790843104117434615624121937495622372_qp &
             )))))))))))))))
 elseif(m.le.0.6_qp) then
         mx=0.45_qp-mc
         elb=     0.8554696151564199914087224774321783838373_qp &
             +mx*(0.1708960726897395844132234165994754905373_qp &
             +mx*(0.1213352290269482260207667564010437464156_qp &
             +mx*(0.1282018835749474096272901529341076494573_qp &
             +mx*(0.1646872814515275597348427294090563472179_qp &
             +mx*(0.2374189087493817423375114793658754489958_qp &
             +mx*(0.3692081047164954516884561039890508294508_qp &
             +mx*(0.6056587338479277173311618664015401963868_qp &
             +mx*(1.0337055615578127436826717513776452476106_qp &
             +mx*(1.8189884893632678849599091011718520567105_qp &
             +mx*(3.2793776512738509375806561547016925831128_qp &
             +mx*(6.0298883807175363312261449542978750456611_qp &
             +mx*(11.269796855577941715109155203721740735793_qp &
             +mx*(21.354577850382834496786315532111529462693_qp &
             )))))))))))))
         eld=     1.04345529511513353426326823569160142342838_qp &
             +mx*(0.77962572192850485048535711388072271372632_qp &
             +mx*(1.02974236093206758187389128668777397528702_qp &
             +mx*(1.62203722341135313022433907993860147395972_qp &
             +mx*(2.78798953118534762046989770119382209443756_qp &
             +mx*(5.04838148737206914685643655935236541332892_qp &
             +mx*(9.46327761194348429539987572314952029503864_qp &
             +mx*(18.1814899494276679043749394081463811247757_qp &
             +mx*(35.5809805911791687037085198750213045708148_qp &
             +mx*(70.6339354619144501276254906239838074917358_qp &
             +mx*(141.828580083433059297030133195739832297859_qp &
             +mx*(287.448751250132166257642182637978103762677_qp &
             +mx*(587.115384649923076181773192202238389711345_qp &
             +mx*(1207.06543522548061603657141890778290399603_qp &
             +mx*(2495.58872724866422273012188618178997342537_qp &
             +mx*(5184.69242939480644062471334944523925163600_qp &
             +mx*(10817.2133369041327524988910635205356016939_qp &
             ))))))))))))))))
 elseif(m.le.0.7_qp) then
         mx=0.35_qp-mc
         elb=     0.8739200618486431359820482173294324246058_qp &
             +mx*(0.1998140574823769459497418213885348159654_qp &
             +mx*(0.1727696158780152128147094051876565603862_qp &
             +mx*(0.2281069132842021671319791750725846795701_qp &
             +mx*(0.3704681411180712197627619157146806221767_qp &
             +mx*(0.6792712528848205545443855883980014994450_qp &
             +mx*(1.3480084966817573020596179874311042267679_qp &
             +mx*(2.8276709768538207038046918250872679902352_qp &
             +mx*(6.1794682501239140840906583219887062092430_qp &
             +mx*(13.935686010342811497608625663457407447757_qp &
             +mx*(32.218929281059722026322932181420383764028_qp &
             +mx*(76.006962959226101026399085304912635262362_qp &
             +mx*(182.32144908775406957609058046006949657416_qp &
             +mx*(443.51507644112648158679360783118806161062_qp &
             +mx*(1091.8547229028388292980623647414961662223_qp &
             +mx*(2715.7658664038195881056269799613407111521_qp &
             )))))))))))))))
         eld=     1.13367833657573316571671258513452768536080_qp &
             +mx*(1.04864317372997039116746991765351150490010_qp &
             +mx*(1.75346504119846451588826580872136305225406_qp &
             +mx*(3.52318272680338551269021618722443199230946_qp &
             +mx*(7.74947641381397458240336356601913534598302_qp &
             +mx*(17.9864500558507330560532617743406294626849_qp &
             +mx*(43.2559163462326133313977294448984936591235_qp &
             +mx*(106.681534454096017031613223924991564429656_qp &
             +mx*(268.098486573117433951562111736259672695883_qp &
             +mx*(683.624114850289804796762005964155730439745_qp &
             +mx*(1763.49708521918740723028849567007874329637_qp &
             +mx*(4592.37475383116380899419201719007475759114_qp &
             +mx*(12053.4410190488892782190764838488156555734_qp &
             +mx*(31846.6630207420816960681624497373078887317_qp &
             +mx*(84621.2213590568080177035346867495326879117_qp &
             +mx*(225956.423182907889987641304430180593010940_qp &
             +mx*(605941.517281758859958050194535269219533685_qp &
             +mx*(1.63108259953926832083633544697688841456604e6_qp &
             )))))))))))))))))
 elseif(m.le.0.8_qp) then
         mx=0.25_qp-mc
         elb=     0.895902820924731621258525533131864225704_qp &
             +mx*(0.243140003766786661947749288357729051637_qp &
             +mx*(0.273081875594105531575351304277604081620_qp &
             +mx*(0.486280007533573323895498576715458103274_qp &
             +mx*(1.082747437228230914750752674136983406683_qp &
             +mx*(2.743445290986452500459431536349945437824_qp &
             +mx*(7.555817828670234627048618342026400847824_qp &
             +mx*(22.05194082493752427472777448620986154515_qp &
             +mx*(67.15640644740229407624192175802742979626_qp &
             +mx*(211.2722537881770961691291434845898538537_qp &
             +mx*(681.9037843053270682273212958093073895805_qp &
             +mx*(2246.956231592536516768812462150619631201_qp &
             +mx*(7531.483865999711792004783423815426725079_qp &
             +mx*(25608.51260130241579018675054866136922157_qp &
             +mx*(88140.74740089604971425934283371277143256_qp &
             +mx*(306564.4242098446591430938434419151070722_qp &
             +mx*(1.076036077811072193752770590363885180738e6_qp &
             +mx*(3.807218502573632648224286313875985190526e6_qp &
             +mx*(1.356638224422139551020110323739879481042e7_qp &
             ))))))))))))))))))
         eld=     1.26061282657491161418014946566845780315983_qp &
             +mx*(1.54866563808267658056930177790599939977154_qp &
             +mx*(3.55366941187160761540650011660758187283401_qp &
             +mx*(9.90044467610439875577300608183010716301714_qp &
             +mx*(30.3205666174524719862025105898574414438275_qp &
             +mx*(98.1802586588830891484913119780870074464833_qp &
             +mx*(329.771010434557055036273670551546757245808_qp &
             +mx*(1136.65598974289039303581967838947708073239_qp &
             +mx*(3993.83433574622979757935610692842933356144_qp &
             +mx*(14242.7295865552708506232731633468180669284_qp &
             +mx*(51394.7572916887209594591528374806790960057_qp &
             +mx*(187246.702914623152141768788258141932569037_qp &
             +mx*(687653.092375389902708761221294282367947659_qp &
             +mx*(2.54238553565398227033448846432182516906624e6_qp &
             +mx*(9.45378121934749027243313241962076028066811e6_qp &
             +mx*(3.53283630179709170835024033154326126569613e7_qp &
             +mx*(1.32593262383393014923560730485845833322771e8_qp &
             +mx*(4.99544968184054821463279808395426941549833e8_qp &
             +mx*(1.88840934729443872364972817525484292678543e9_qp &
             +mx*(7.16026753447893719179055010636502508063102e9_qp &
             +mx*(2.72233079469633962247554894093591262281929e10_qp &
         ))))))))))))))))))))
 elseif(m.le.0.85_qp) then
         mx=0.175_qp-mc
         elb=     0.915922052601931494319853880201442948834592_qp &
             +mx*(0.294714252429483394379515488141632749820347_qp &
             +mx*(0.435776709264636140422971598963772380161131_qp &
             +mx*(1.067328246493644238508159085364429570207744_qp &
             +mx*(3.327844118563268085074646976514979307993733_qp &
             +mx*(11.90406004445092906188837729711173326621810_qp &
             +mx*(46.47838820224626393512400481776284680677096_qp &
             +mx*(192.7556002578809476962739389101964074608802_qp &
             +mx*(835.3356299261900063712302517586717381557137_qp &
             +mx*(3743.124548343029102644419963712353854902019_qp &
             +mx*(17219.07731004063094108708549153310467326395_qp &
             +mx*(80904.60401669850158353080543152212152282878_qp &
             +mx*(386808.3292751742460123683674607895217760313_qp &
             +mx*(1.876487670110449342170327796786290400635732e6_qp &
             +mx*(9.216559908641567755240142886998737950775910e6_qp &
             ))))))))))))))
         eld=     1.402200569110579095046054435635136986038164_qp &
             +mx*(2.322205897861749446534352741005347103992773_qp &
             +mx*(7.462158366466719682730245467372788273333992_qp &
             +mx*(29.43506890797307903104978364254987042421285_qp &
             +mx*(128.1590924337895775262509354898066132182429_qp &
             +mx*(591.0807036911982326384997979640812493154316_qp &
             +mx*(2830.546229607726377048576057043685514661188_qp &
             +mx*(13917.76431889392229954434840686741305556862_qp &
             +mx*(69786.10525163921228258055074102587429394212_qp &
             +mx*(355234.1420341879634781808899208309503519936_qp &
             +mx*(1.830019186413931053503912913904321703777885e6_qp &
             +mx*(9.519610812032515607466102200648641326190483e6_qp &
             +mx*(4.992086875574849453986274042758566713803723e7_qp &
             +mx*(2.635677009826023473846461512029006874800883e8_qp &
             +mx*(1.399645765120061118824228996253541612110338e9_qp &
             +mx*(7.469935792837635004663183580452618726280406e9_qp &
             +mx*(4.004155595835610574316003488168804738481448e10_qp &
             +mx*(2.154630668144966654449602981243932210695662e11_qp &
             )))))))))))))))))
 else
     mx=0.125_qp-mc
         elb=     0.931906061029524827613331428871579482766771_qp &
             +mx*(0.348448029538453860999386797137074571589376_qp &
             +mx*(0.666809178846938247558793864839434184202736_qp &
             +mx*(2.210769135708128662563678717558631573758222_qp &
             +mx*(9.491765048913406881414290930355300611703187_qp &
             +mx*(47.09304791027740853381457907791343619298913_qp &
             +mx*(255.9200460211233087050940506395442544885608_qp &
             +mx*(1480.029532675805407554800779436693505109703_qp &
             +mx*(8954.040904734313578374783155553041065984547_qp &
             +mx*(56052.48220982686949967604699243627330816542_qp &
             +mx*(360395.7241626000916973524840479780937869149_qp &
             +mx*(2.367539415273216077520928806581689330885103e6_qp &
             +mx*(1.582994957277684102454906900025484391190264e7_qp &
             +mx*(1.074158093278511100137056972128875270067228e8_qp &
             +mx*(7.380585460239595691878086073095523043390649e8_qp &
             +mx*(5.126022002555101496684687154904781856830296e9_qp &
             +mx*(3.593534065502416588712409180013118409428367e10_qp &
             +mx*(2.539881257612812212004146637239987308133582e11_qp &
             +mx*(1.808180007145359569674767150594344316702507e12_qp &
         ))))))))))))))))))
         eld=     1.541690112721819084362258323861459983048179_qp &
             +mx*(3.379176214579645449453938918349243359477706_qp &
             +mx*(14.94058385670236671625328259137998668324435_qp &
             +mx*(81.91773929235074880784578753539752529822986_qp &
             +mx*(497.4900546551479866036061853049402721939835_qp &
             +mx*(3205.184010234846235275447901572262470252768_qp &
             +mx*(21457.32237355321925571253220641357074594515_qp &
             +mx*(147557.0156564174712105689758692929775004292_qp &
             +mx*(1.035045290185256525452269053775273002725343e6_qp &
             +mx*(7.371922334832212125197513363695905834126154e6_qp &
             +mx*(5.314344395142401141792228169170505958906345e7_qp &
             +mx*(3.868823475795976312985118115567305767603128e8_qp &
             +mx*(2.839458401528033778425531336599562337200510e9_qp &
             +mx*(2.098266122943898941547136470383199468548861e10_qp &
             +mx*(1.559617754017662417944194874282275405676282e11_qp &
             +mx*(1.165096220419884791236699872205721392201682e12_qp &
             +mx*(8.742012983013913804987431275193291316808766e12_qp &
             +mx*(6.584725462672366918676967847406180155459650e13_qp &
             +mx*(4.976798737062434393396993620379481464465749e14_qp &
             +mx*(3.773018634056605404718444239040628892506293e15_qp &
             +mx*(2.868263194837819660109735981973458220407767e16_qp &
         ))))))))))))))))))))
 endif
 mcold=mc
 elbold=elb
 eldold=eld
 return
 end

 elemental subroutine serbd_qp(y,m,b,d)
 !
 ! Simultaneous computation of associate elliptic integrals,
 ! B(phi|m) and D(phi|m), for small arguments by the series expansion
 !
 ! Reference: Fukushima, T (2012) J. Comp. Appl. Math., 235, 4140-4148
 !   Precise and fast computation of general incomplete elliptic integral
 !   of second kind by half and double argument transformations
 !
 real(qp),intent(in) :: y,m !!SJT
 real(qp),intent(out) :: b,d !!SJT
 !
 !real(qp) F1,F2,F3,F4 !!SJT: OG
 real(qp) F1,F2,F3,F4,F5,F6,F7,F8,F9,FA,FB !!SJT: non-parameter variables
 !real(qp) F10,F20,F21,F30,F31,F40,F41,F42
 !real(qp) F5,F50,F51,F52,F6,F60,F61,F62,F63
 !real(qp) F7,F70,F71,F72,F73,F8,F80,F81,F82,F83,F84
 !real(qp) F9,F90,F91,F92,F93,F94
 !real(qp) FA,FA0,FA1,FA2,FA3,FA4,FA5
 !real(qp) FB,FB0,FB1,FB2,FB3,FB4,FB5
 real(qp),parameter :: F10=1._qp/6._qp
 real(qp),parameter :: F20=3._qp/40._qp
 real(qp),parameter :: F21=2._qp/40._qp
 real(qp),parameter :: F30=5._qp/112._qp
 real(qp),parameter :: F31=3._qp/112._qp
 real(qp),parameter :: F40=35._qp/1152._qp
 real(qp),parameter :: F41=20._qp/1152._qp
 real(qp),parameter :: F42=18._qp/1152._qp
 real(qp),parameter :: F50=63._qp/2816._qp
 real(qp),parameter :: F51=35._qp/2816._qp
 real(qp),parameter :: F52=30._qp/2816._qp
 real(qp),parameter :: F60=231._qp/13312._qp
 real(qp),parameter :: F61=126._qp/13312._qp
 real(qp),parameter :: F62=105._qp/13312._qp
 real(qp),parameter :: F63=100._qp/13312._qp
 real(qp),parameter :: F70=429._qp/30720._qp
 real(qp),parameter :: F71=231._qp/30720._qp
 real(qp),parameter :: F72=189._qp/30720._qp
 real(qp),parameter :: F73=175._qp/30720._qp
 real(qp),parameter :: F80=6435._qp/557056._qp
 real(qp),parameter :: F81=3432._qp/557056._qp
 real(qp),parameter :: F82=2722._qp/557056._qp
 real(qp),parameter :: F83=2520._qp/557056._qp
 real(qp),parameter :: F84=2450._qp/557056._qp
 real(qp),parameter :: F90=12155._qp/1245184._qp
 real(qp),parameter :: F91=6435._qp/1245184._qp
 real(qp),parameter :: F92=5148._qp/1245184._qp
 real(qp),parameter :: F93=4620._qp/1245184._qp
 real(qp),parameter :: F94=4410._qp/1245184._qp
 real(qp),parameter :: FA0=46189._qp/5505024._qp
 real(qp),parameter :: FA1=24310._qp/5505024._qp
 real(qp),parameter :: FA2=19305._qp/5505024._qp
 real(qp),parameter :: FA3=17160._qp/5505024._qp
 real(qp),parameter :: FA4=16170._qp/5505024._qp
 real(qp),parameter :: FA5=15876._qp/5505024._qp
 real(qp),parameter :: FB0=88179._qp/12058624._qp
 real(qp),parameter :: FB1=46189._qp/12058624._qp
 real(qp),parameter :: FB2=36465._qp/12058624._qp
 real(qp),parameter :: FB3=32175._qp/12058624._qp
 real(qp),parameter :: FB4=30030._qp/12058624._qp
 real(qp),parameter :: FB5=29106._qp/12058624._qp
 !real(qp) A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB
 real(qp),parameter :: A1=3._qp/5._qp
 real(qp),parameter :: A2=5._qp/7._qp
 real(qp),parameter :: A3=7._qp/9._qp
 real(qp),parameter :: A4=9._qp/11._qp
 real(qp),parameter :: A5=11._qp/13._qp
 real(qp),parameter :: A6=13._qp/15._qp
 real(qp),parameter :: A7=15._qp/17._qp
 real(qp),parameter :: A8=17._qp/19._qp
 real(qp),parameter :: A9=19._qp/21._qp
 real(qp),parameter :: AA=21._qp/23._qp
 real(qp),parameter :: AB=23._qp/25._qp
 real(qp) B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB
 !real(qp) D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,DA,DB !!SJT: OG
 real(qp) D1,D2,D3,D4,D5,D6,D7,D8,D9,DA,DB !!SJT: non-parameter variables
 real(qp),parameter :: D0=1._qp/3._qp
 !
 !        write(*,*) "(serbd) y,m=",y,m
 F1=F10+m*F10
 F2=F20+m*(F21+m*F20)
 F3=F30+m*(F31+m*(F31+m*F30))
 F4=F40+m*(F41+m*(F42+m*(F41+m*F40)))
 F5=F50+m*(F51+m*(F52+m*(F52+m*(F51+m*F50))))
 F6=F60+m*(F61+m*(F62+m*(F63+m*(F62+m*(F61+m*F60)))))
 F7=F70+m*(F71+m*(F72+m*(F73+m*(F73+m*(F72+m*(F71+m*F70))))))
 F8=F80+m*(F81+m*(F82+m*(F83+m*(F84+m*(F83+m*(F82+m*(F81+m*F80)))))))
 F9=F90+m*(F91+m*(F92+m*(F93+m*(F94+m*(F94+m*(F93 &
     +m*(F92+m*(F91+m*F90))))))))
 FA=FA0+m*(FA1+m*(FA2+m*(FA3+m*(FA4+m*(FA5+m*(FA4 &
     +m*(FA3+m*(FA2+m*(FA1+m*FA0)))))))))
 FB=FB0+m*(FB1+m*(FB2+m*(FB3+m*(FB4+m*(FB5+m*(FB5+ &
     m*(FB4+m*(FB3+m*(FB2+m*(FB1+m*FB0))))))))))
 !
 D1=F1*A1
 D2=F2*A2
 D3=F3*A3
 D4=F4*A4
 D5=F5*A5
 D6=F6*A6
 D7=F7*A7
 D8=F8*A8
 D9=F9*A9
 DA=FA*AA
 DB=FB*AB
 d=D0+y*(D1+y*(D2+y*(D3+y*(D4+y*(D5+y*(D6+y*(D7+y*(D8 &
     +y*(D9+y*(DA+y*DB))))))))))
 B1=F1-D0
 B2=F2-D1
 B3=F3-D2
 B4=F4-D3
 B5=F5-D4
 B6=F6-D5
 B7=F7-D6
 B8=F8-D7
 B9=F9-D8
 BA=FA-D9
 BB=FB-DA
 b=1._qp+y*(B1+y*(B2+y*(B3+y*(B4+y*(B5+y*(B6+y*(B7+y*(B8 &
     +y*(B9+y*(BA+y*BB))))))))))
 return
 end

 elemental real(qp) function serj_qp(y,n,m) result(serj)
 !
 ! Computation of associate elliptic integral J(phi,n|m)
 ! for small arguments by the series expansion
 !
 ! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !   Precise and fast computation of a general incomplete elliptic integral
 !   of third kind by half and double argument transformations
 !
 real(qp),intent(in) :: y,n,m !!SJT
 
 real(qp) J1,J2,J3,J4,J5,J6,J7,J8,J9,JA
 
 !real(qp) J100,J200,J201,J210,J300,J301,J302,J310,J311,J320
 !real(qp) J400,J401,J402,J403,J410,J411,J412,J420,J421,J430
 !real(qp) J500,J501,J502,J503,J504,J510,J511,J512,J513,J520
 !real(qp) J521,J522,J530,J531,J540
 !real(qp) J600,J601,J602,J603,J604,J605,J610,J611,J612,J613,J614
 !real(qp) J620,J621,J622,J623,J630,J631,J632,J640,J641,J650
 !real(qp) J700,J701,J702,J703,J704,J705,J706
 !real(qp) J710,J711,J712,J713,J714,J715,J720,J721,J722,J723,J724
 !real(qp) J730,J731,J732,J733,J740,J741,J742,J750,J751,J760
 !real(qp) J800,J801,J802,J803,J804,J805,J806,J807
 !real(qp) J810,J811,J812,J813,J814,J815,J816
 !real(qp) J820,J821,J822,J823,J824,J825,J830,J831,J832,J833,J834
 !real(qp) J840,J841,J842,J843,J850,J851,J852,J860,J861,J870
 !real(qp) J900,J901,J902,J903,J904,J905,J906,J907,J908
 !real(qp) J910,J911,J912,J913,J914,J915,J916,J917
 !real(qp) J920,J921,J922,J923,J924,J925,J926
 !real(qp) J930,J931,J932,J933,J934,J935,J940,J941,J942,J943,J944
 !real(qp) J950,J951,J952,J953,J960,J961,J962,J970,J971,J980
 !real(qp) JA00,JA01,JA02,JA03,JA04,JA05,JA06,JA07,JA08,JA09
 !real(qp) JA10,JA11,JA12,JA13,JA14,JA15,JA16,JA17,JA18
 !real(qp) JA20,JA21,JA22,JA23,JA24,JA25,JA26,JA27
 !real(qp) JA30,JA31,JA32,JA33,JA34,JA35,JA36
 !real(qp) JA40,JA41,JA42,JA43,JA44,JA45,JA50,JA51,JA52,JA53,JA54
 !real(qp) JA60,JA61,JA62,JA63,JA70,JA71,JA72,JA80,JA81,JA90
 
 real(qp),parameter :: J100=1._qp/3._qp
 
 real(qp),parameter :: J200=1._qp/10._qp
 real(qp),parameter :: J201=2._qp/10._qp
 real(qp),parameter :: J210=1._qp/10._qp
 
 real(qp),parameter :: J300=3._qp/56._qp
 real(qp),parameter :: J301=4._qp/56._qp
 real(qp),parameter :: J302=8._qp/56._qp
 real(qp),parameter :: J310=2._qp/56._qp
 real(qp),parameter :: J311=4._qp/56._qp
 real(qp),parameter :: J320=3._qp/56._qp
 
 real(qp),parameter :: J400=5._qp/144._qp
 real(qp),parameter :: J401=6._qp/144._qp
 real(qp),parameter :: J402=8._qp/144._qp
 real(qp),parameter :: J403=16._qp/144._qp
 real(qp),parameter :: J410=3._qp/144._qp
 real(qp),parameter :: J411=4._qp/144._qp
 real(qp),parameter :: J412=8._qp/144._qp
 real(qp),parameter :: J420=3._qp/144._qp
 real(qp),parameter :: J421=6._qp/144._qp
 real(qp),parameter :: J430=5._qp/144._qp
 
 real(qp),parameter :: J500=35._qp/1408._qp
 real(qp),parameter :: J501=40._qp/1408._qp
 real(qp),parameter :: J502=48._qp/1408._qp
 real(qp),parameter :: J503=64._qp/1408._qp
 real(qp),parameter :: J504=128._qp/1408._qp
 real(qp),parameter :: J510=20._qp/1408._qp
 real(qp),parameter :: J511=24._qp/1408._qp
 real(qp),parameter :: J512=32._qp/1408._qp
 real(qp),parameter :: J513=64._qp/1408._qp
 real(qp),parameter :: J520=18._qp/1408._qp
 real(qp),parameter :: J521=24._qp/1408._qp
 real(qp),parameter :: J522=48._qp/1408._qp
 real(qp),parameter :: J530=20._qp/1408._qp
 real(qp),parameter :: J531=40._qp/1408._qp
 real(qp),parameter :: J540=35._qp/1408._qp
 
 real(qp),parameter :: J600=63._qp/3328._qp
 real(qp),parameter :: J601=70._qp/3328._qp
 real(qp),parameter :: J602=80._qp/3328._qp
 real(qp),parameter :: J603=96._qp/3328._qp
 real(qp),parameter :: J604=128._qp/3328._qp
 real(qp),parameter :: J605=256._qp/3328._qp
 real(qp),parameter :: J610=35._qp/3328._qp
 real(qp),parameter :: J611=40._qp/3328._qp
 real(qp),parameter :: J612=48._qp/3328._qp
 real(qp),parameter :: J613=64._qp/3328._qp
 real(qp),parameter :: J614=128._qp/3328._qp
 real(qp),parameter :: J620=30._qp/3328._qp
 real(qp),parameter :: J621=36._qp/3328._qp
 real(qp),parameter :: J622=48._qp/3328._qp
 real(qp),parameter :: J623=96._qp/3328._qp
 real(qp),parameter :: J630=30._qp/3328._qp
 real(qp),parameter :: J631=40._qp/3328._qp
 real(qp),parameter :: J632=80._qp/3328._qp
 real(qp),parameter :: J640=35._qp/3328._qp
 real(qp),parameter :: J641=70._qp/3328._qp
 real(qp),parameter :: J650=63._qp/3328._qp
 
 real(qp),parameter :: J700=231._qp/15360._qp
 real(qp),parameter :: J701=252._qp/15360._qp
 real(qp),parameter :: J702=280._qp/15360._qp
 real(qp),parameter :: J703=320._qp/15360._qp
 real(qp),parameter :: J704=384._qp/15360._qp
 real(qp),parameter :: J705=512._qp/15360._qp
 real(qp),parameter :: J706=1024._qp/15360._qp
 real(qp),parameter :: J710=126._qp/15360._qp
 real(qp),parameter :: J711=140._qp/15360._qp
 real(qp),parameter :: J712=160._qp/15360._qp
 real(qp),parameter :: J713=192._qp/15360._qp
 real(qp),parameter :: J714=256._qp/15360._qp
 real(qp),parameter :: J715=512._qp/15360._qp
 real(qp),parameter :: J720=105._qp/15360._qp
 real(qp),parameter :: J721=120._qp/15360._qp
 real(qp),parameter :: J722=144._qp/15360._qp
 real(qp),parameter :: J723=192._qp/15360._qp
 real(qp),parameter :: J724=384._qp/15360._qp
 real(qp),parameter :: J730=100._qp/15360._qp
 real(qp),parameter :: J731=120._qp/15360._qp
 real(qp),parameter :: J732=160._qp/15360._qp
 real(qp),parameter :: J733=320._qp/15360._qp
 real(qp),parameter :: J740=105._qp/15360._qp
 real(qp),parameter :: J741=140._qp/15360._qp
 real(qp),parameter :: J742=280._qp/15360._qp
 real(qp),parameter :: J750=126._qp/15360._qp
 real(qp),parameter :: J751=252._qp/15360._qp
 real(qp),parameter :: J760=231._qp/15360._qp
 
 real(qp),parameter :: J800=429._qp/34816._qp
 real(qp),parameter :: J801=462._qp/34816._qp
 real(qp),parameter :: J802=504._qp/34816._qp
 real(qp),parameter :: J803=560._qp/34816._qp
 real(qp),parameter :: J804=640._qp/34816._qp
 real(qp),parameter :: J805=768._qp/34816._qp
 real(qp),parameter :: J806=1024._qp/34816._qp
 real(qp),parameter :: J807=2048._qp/34816._qp
 real(qp),parameter :: J810=231._qp/34816._qp
 real(qp),parameter :: J811=252._qp/34816._qp
 real(qp),parameter :: J812=280._qp/34816._qp
 real(qp),parameter :: J813=320._qp/34816._qp
 real(qp),parameter :: J814=384._qp/34816._qp
 real(qp),parameter :: J815=512._qp/34816._qp
 real(qp),parameter :: J816=1024._qp/34816._qp
 real(qp),parameter :: J820=189._qp/34816._qp
 real(qp),parameter :: J821=210._qp/34816._qp
 real(qp),parameter :: J822=240._qp/34816._qp
 real(qp),parameter :: J823=288._qp/34816._qp
 real(qp),parameter :: J824=284._qp/34816._qp
 real(qp),parameter :: J825=768._qp/34816._qp
 real(qp),parameter :: J830=175._qp/34816._qp
 real(qp),parameter :: J831=200._qp/34816._qp
 real(qp),parameter :: J832=240._qp/34816._qp
 real(qp),parameter :: J833=320._qp/34816._qp
 real(qp),parameter :: J834=640._qp/34816._qp
 real(qp),parameter :: J840=175._qp/34816._qp
 real(qp),parameter :: J841=210._qp/34816._qp
 real(qp),parameter :: J842=280._qp/34816._qp
 real(qp),parameter :: J843=560._qp/34816._qp
 real(qp),parameter :: J850=189._qp/34816._qp
 real(qp),parameter :: J851=252._qp/34816._qp
 real(qp),parameter :: J852=504._qp/34816._qp
 real(qp),parameter :: J860=231._qp/34816._qp
 real(qp),parameter :: J861=462._qp/34816._qp
 real(qp),parameter :: J870=429._qp/34816._qp
 
 real(qp),parameter :: J900=6435._qp/622592._qp
 real(qp),parameter :: J901=6864._qp/622592._qp
 real(qp),parameter :: J902=7392._qp/622592._qp
 real(qp),parameter :: J903=8064._qp/622592._qp
 real(qp),parameter :: J904=8960._qp/622592._qp
 real(qp),parameter :: J905=10240._qp/622592._qp
 real(qp),parameter :: J906=12288._qp/622592._qp
 real(qp),parameter :: J907=16384._qp/622592._qp
 real(qp),parameter :: J908=32768._qp/622592._qp
 real(qp),parameter :: J910=3432._qp/622592._qp
 real(qp),parameter :: J911=3696._qp/622592._qp
 real(qp),parameter :: J912=4032._qp/622592._qp
 real(qp),parameter :: J913=4480._qp/622592._qp
 real(qp),parameter :: J914=5120._qp/622592._qp
 real(qp),parameter :: J915=6144._qp/622592._qp
 real(qp),parameter :: J916=8192._qp/622592._qp
 real(qp),parameter :: J917=16384._qp/622592._qp
 real(qp),parameter :: J920=2772._qp/622592._qp
 real(qp),parameter :: J921=3024._qp/622592._qp
 real(qp),parameter :: J922=3360._qp/622592._qp
 real(qp),parameter :: J923=3840._qp/622592._qp
 real(qp),parameter :: J924=4608._qp/622592._qp
 real(qp),parameter :: J925=6144._qp/622592._qp
 real(qp),parameter :: J926=12288._qp/622592._qp
 real(qp),parameter :: J930=2520._qp/622592._qp
 real(qp),parameter :: J931=2800._qp/622592._qp
 real(qp),parameter :: J932=3200._qp/622592._qp
 real(qp),parameter :: J933=3840._qp/622592._qp
 real(qp),parameter :: J934=5120._qp/622592._qp
 real(qp),parameter :: J935=10240._qp/622592._qp
 real(qp),parameter :: J940=2450._qp/622592._qp
 real(qp),parameter :: J941=2800._qp/622592._qp
 real(qp),parameter :: J942=3360._qp/622592._qp
 real(qp),parameter :: J943=4480._qp/622592._qp
 real(qp),parameter :: J944=8960._qp/622592._qp
 real(qp),parameter :: J950=2520._qp/622592._qp
 real(qp),parameter :: J951=3024._qp/622592._qp
 real(qp),parameter :: J952=4032._qp/622592._qp
 real(qp),parameter :: J953=8064._qp/622592._qp
 real(qp),parameter :: J960=2772._qp/622592._qp
 real(qp),parameter :: J961=3696._qp/622592._qp
 real(qp),parameter :: J962=7392._qp/622592._qp
 real(qp),parameter :: J970=3432._qp/622592._qp
 real(qp),parameter :: J971=6864._qp/622592._qp
 real(qp),parameter :: J980=6435._qp/622592._qp
 
 real(qp),parameter :: JA00=12155._qp/1376256._qp
 real(qp),parameter :: JA01=12870._qp/1376256._qp
 real(qp),parameter :: JA02=13728._qp/1376256._qp
 real(qp),parameter :: JA03=14784._qp/1376256._qp
 real(qp),parameter :: JA04=16128._qp/1376256._qp
 real(qp),parameter :: JA05=17920._qp/1376256._qp
 real(qp),parameter :: JA06=20480._qp/1376256._qp
 real(qp),parameter :: JA07=24576._qp/1376256._qp
 real(qp),parameter :: JA08=32768._qp/1376256._qp
 real(qp),parameter :: JA09=65536._qp/1376256._qp
 real(qp),parameter :: JA10=6435._qp/1376256._qp
 real(qp),parameter :: JA11=6864._qp/1376256._qp
 real(qp),parameter :: JA12=7392._qp/1376256._qp
 real(qp),parameter :: JA13=8064._qp/1376256._qp
 real(qp),parameter :: JA14=8960._qp/1376256._qp
 real(qp),parameter :: JA15=10240._qp/1376256._qp
 real(qp),parameter :: JA16=12288._qp/1376256._qp
 real(qp),parameter :: JA17=16384._qp/1376256._qp
 real(qp),parameter :: JA18=32768._qp/1376256._qp
 real(qp),parameter :: JA20=5148._qp/1376256._qp
 real(qp),parameter :: JA21=5544._qp/1376256._qp
 real(qp),parameter :: JA22=6048._qp/1376256._qp
 real(qp),parameter :: JA23=6720._qp/1376256._qp
 real(qp),parameter :: JA24=7680._qp/1376256._qp
 real(qp),parameter :: JA25=9216._qp/1376256._qp
 real(qp),parameter :: JA26=12288._qp/1376256._qp
 real(qp),parameter :: JA27=24576._qp/1376256._qp
 real(qp),parameter :: JA30=4620._qp/1376256._qp
 real(qp),parameter :: JA31=5040._qp/1376256._qp
 real(qp),parameter :: JA32=5600._qp/1376256._qp
 real(qp),parameter :: JA33=6400._qp/1376256._qp
 real(qp),parameter :: JA34=7680._qp/1376256._qp
 real(qp),parameter :: JA35=10240._qp/1376256._qp
 real(qp),parameter :: JA36=20480._qp/1376256._qp
 real(qp),parameter :: JA40=4410._qp/1376256._qp
 real(qp),parameter :: JA41=4900._qp/1376256._qp
 real(qp),parameter :: JA42=5600._qp/1376256._qp
 real(qp),parameter :: JA43=6720._qp/1376256._qp
 real(qp),parameter :: JA44=8960._qp/1376256._qp
 real(qp),parameter :: JA45=17920._qp/1376256._qp
 real(qp),parameter :: JA50=4410._qp/1376256._qp
 real(qp),parameter :: JA51=5040._qp/1376256._qp
 real(qp),parameter :: JA52=6048._qp/1376256._qp
 real(qp),parameter :: JA53=8064._qp/1376256._qp
 real(qp),parameter :: JA54=16128._qp/1376256._qp
 real(qp),parameter :: JA60=4620._qp/1376256._qp
 real(qp),parameter :: JA61=5544._qp/1376256._qp
 real(qp),parameter :: JA62=7392._qp/1376256._qp
 real(qp),parameter :: JA63=14784._qp/1376256._qp
 real(qp),parameter :: JA70=5148._qp/1376256._qp
 real(qp),parameter :: JA71=6864._qp/1376256._qp
 real(qp),parameter :: JA72=13728._qp/1376256._qp
 real(qp),parameter :: JA80=6435._qp/1376256._qp
 real(qp),parameter :: JA81=12870._qp/1376256._qp
 real(qp),parameter :: JA90=12155._qp/1376256._qp
 
 ! write(*,"(a20,1p3e10.2)") "(serj) y,n,m=",y,n,m
 
 J1=J100
 J2=J200+n*J201+m*J210
 !if(y.lt.9.2925066E-09) then
 !    serj=y*(J1+y*J2)
 ! write(*,"(a20,1pe10.2)") "(serj) J2",serj
 !    return
 !endif
 
 J3=J300+n*(J301+n*J302)+m*(J310+n*J311+m*J320)
 !if(y.lt.4.3667383d-06) then
 !    serj=y*(J1+y*(J2+y*J3))
 ! write(*,"(a20,1pe10.2)") "(serj) J3",serj
 !    return
 !endif
 
 J4=J400+n*(J401+n*(J402+n*J403)) &
     +m*(J410+n*(J411+n*J412)+m*(J420+n*J421+m*J430))
 !if(y.le.9.4990006E-05) then
 !    serj=y*(J1+y*(J2+y*(J3+y*J4)))
 ! write(*,"(a20,1pe10.2)") "(serj) J4",serj
 !    return
 !endif
 
 J5=J500+n*(J501+n*(J502+n*(J503+n*J504))) &
     +m*(J510+n*(J511+n*(J512+n*J513)) &
     +m*(J520+n*(J521+n*J522)+m*(J530+n*J531+m*J540)))
 if(y.le.6.0369310e-04_qp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*J5))))
 ! write(*,"(a20,1pe10.2)") "(serj) J5",serj
     return
 endif
 
 J6=J600+n*(J601+n*(J602+n*(J603+n*(J604+n*J605)))) &
     +m*(J610+n*(J611+n*(J612+n*(J613+n*J614))) &
     +m*(J620+n*(J621+n*(J622+n*J623)) &
     +m*(J630+n*(J631+n*J632)+m*(J640+n*J641+m*J650))))
 if(y.le.2.0727505e-03_qp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*J6)))))
 ! write(*,"(a20,1pe10.2)") "(serj) J6",serj
     return
 endif
 
 J7=J700+n*(J701+n*(J702+n*(J703+n*(J704+n*(J705+n*J706))))) &
     +m*(J710+n*(J711+n*(J712+n*(J713+n*(J714+n*J715)))) &
     +m*(J720+n*(J721+n*(J722+n*(J723+n*J724))) &
     +m*(J730+n*(J731+n*(J732+n*J733)) &
     +m*(J740+n*(J741+n*J742)+m*(J750+n*J751+m*J760)))))
 if(y.le.5.0047026e-03_qp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*J7))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J7",serj
    return
 endif
 
 J8=J800+n*(J801+n*(J802+n*(J803+n*(J804+n*(J805+n*(J806+n*J807)))))) &
     +m*(J810+n*(J811+n*(J812+n*(J813+n*(J814+n*(J815+n*J816))))) &
     +m*(J820+n*(J821+n*(J822+n*(J823+n*(J824+n*J825)))) &
     +m*(J830+n*(J831+n*(J832+n*(J833+n*J834))) &
     +m*(J840+n*(J841+n*(J842+n*J843)) &
     +m*(J850+n*(J851+n*J852)+m*(J860+n*J861+m*J870))))))
 if(y.le.9.6961652e-03_qp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*J8)))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J8",serj
     return
 endif
 
 J9=J900+n*(J901+n*(J902+n*(J903+n*(J904+n*(J905+n*(J906+n*(J907+n*J908))))))) &
     +m*(J910+n*(J911+n*(J912+n*(J913+n*(J914+n*(J915+n*(J916+n*J917)))))) &
     +m*(J920+n*(J921+n*(J922+n*(J923+n*(J924+n*(J925+n*J926))))) &
     +m*(J930+n*(J931+n*(J932+n*(J933+n*(J934+n*J935)))) &
     +m*(J940+n*(J941+n*(J942+n*(J943+n*J944))) &
     +m*(J950+n*(J951+n*(J952+n*J953)) &
     +m*(J960+n*(J961+n*J962)+m*(J970+n*J971+m*J980)))))))
 if(y.le.1.6220210e-02_qp) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*(J8+y*J9))))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J9",serj
     return
 endif
 
 JA=JA00+n*(JA01+n*(JA02+n*(JA03+n*(JA04+n*(JA05+n*(JA06+n*(JA07+n*(JA08+n*JA09)))))))) &
     +m*(JA10+n*(JA11+n*(JA12+n*(JA13+n*(JA14+n*(JA15+n*(JA16+n*(JA17+n*JA18))))))) &
     +m*(JA20+n*(JA21+n*(JA22+n*(JA23+n*(JA24+n*(JA25+n*(JA26+n*JA27)))))) &
     +m*(JA30+n*(JA31+n*(JA32+n*(JA33+n*(JA34+n*(JA35+n*JA36))))) &
     +m*(JA40+n*(JA41+n*(JA42+n*(JA43+n*(JA44+n*JA45)))) &
     +m*(JA50+n*(JA51+n*(JA52+n*(JA53+n*JA54))) &
     +m*(JA60+n*(JA61+n*(JA62+n*JA63)) &
     +m*(JA70+n*(JA71+n*JA72)+m*(JA80+n*JA81+m*JA90))))))))
 !if(y.lt.2.4482875d-02) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*(J8+y*(J9+y*JA)))))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J10",serj
     return
 !endif
 
 end

 elemental real(qp) function uatan_qp(t,h) result(uatan)
 !
 ! Universal arctangent function
 !
 ! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !   Precise and fast computation of a general incomplete elliptic integral
 !   of third kind by half and double argument transformations
 !
 real(qp),intent(in) :: t,h !!SJT
 real(qp) :: z,y,x !!SJT
 real(qp) a,r,ri,hold,rold,riold
 !real(qp) A3,A5,A7,A9,A11,A13,A15,A17,A19,A21,A23,A25
 !data hold/1.d0/, rold/1.d0/,riold/1.d0/ !!SJT: disable for thread safety
 !save hold,rold,riold !!SJT: disable for thread safety
 real(qp),parameter :: A3=1._qp/3._qp
 real(qp),parameter :: A5=1._qp/5._qp
 real(qp),parameter :: A7=1._qp/7._qp
 real(qp),parameter :: A9=1._qp/9._qp
 real(qp),parameter :: A11=1._qp/11._qp
 real(qp),parameter :: A13=1._qp/13._qp
 real(qp),parameter :: A15=1._qp/15._qp
 real(qp),parameter :: A17=1._qp/17._qp
 real(qp),parameter :: A19=1._qp/19._qp
 real(qp),parameter :: A21=1._qp/21._qp
 real(qp),parameter :: A23=1._qp/23._qp
 real(qp),parameter :: A25=1._qp/25._qp
 !!!!!SJT: replace data statment with manual value assignment for thread safety
 hold=1._qp; rold=1._qp; riold=1._qp
 !!!!!SJT: end replace data statment with manual value assignment for thread safety
 !
 ! write(*,*) "(uatan) t,h=",t,h
 !
 z=-h*t*t
 a=abs(z)
 !
 ! write(*,*) "(uatan) z=",z
 if(a.lt.3.3306691e-16_qp) then
     uatan=t
 ! write(*,*) "0"
 elseif(a.lt.2.3560805e-08_qp) then
     uatan=t*(1.e0_qp+z*A3)
 ! write(*,*) "1"
 elseif(a.lt.9.1939631e-06_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*A5))
 ! write(*,*) "2"
 elseif(a.lt.1.7779240e-04_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*A7)))
 ! write(*,*) "3"
 elseif(a.lt.1.0407839e-03_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*(A7+z*A9))))
 ! write(*,*) "4"
 elseif(a.lt.3.3616998e-03_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*(A7+z*(A9+z*A11)))))
 ! write(*,*) "5"
 elseif(a.lt.7.7408014e-03_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*A13))))))
 ! write(*,*) "6"
 elseif(a.lt.1.4437181e-02_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*A15)))))))
 ! write(*,*) "7"
 elseif(a.lt.2.3407312e-02_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*A17))))))))
 ! write(*,*) "8"
 elseif(a.lt.3.4416203e-02_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*A19)))))))))
 ! write(*,*) "9"
 elseif(z.lt.0.e0_qp) then
     if(abs(h-hold).lt.1.e-16_qp) then
         r=rold; ri=riold
     else
         r=sqrt(h); ri=1.e0_qp/r; hold=h; rold=r; riold=ri
     endif
     uatan=atan(r*t)*ri
 ! write(*,*) "atan"
 elseif(a.lt.4.7138547e-02_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*A21))))))))))
 ! write(*,*) "A"
 elseif(a.lt.6.1227405e-02_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*(A21+z*A23)))))))))))
 ! write(*,*) "B"
 elseif(a.lt.7.6353468e-02_qp) then
     uatan=t*(1.e0_qp+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*(A21+z*(A23+z*A25))))))))))))
 ! write(*,*) "C"
 else
     if(abs(h-hold).lt.1.e-16_qp) then
         r=rold; ri=riold
     else
         r=sqrt(-h); ri=1.e0_qp/r; hold=h; rold=r; riold=ri
     endif
     y=r*t
     x=log((1.e0_qp+y)/(1.e0_qp-y))*0.5e0_qp
     uatan=x*ri
 ! write(*,*) "hyper"
 endif
 !
 ! write(*,*) "(uatan) uatan=",uatan
 !
 return
 end
end module xelbdj2_all_routines
