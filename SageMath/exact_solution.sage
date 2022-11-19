#Sage script for thermochemical problem: starting from a presumed stream function
###########################################################################BEGIN SETUP
reset() #clear current session
import numpy as np
import matplotlib.pyplot as plt

##Begin Declare Variables
x,z,t=var('x,z,t')
xb=var('xb')
k,zI,zR=var('k,zI,zR')
a,b,c,d=var('a,b,c,d')
RaT,RaC=var('RaT,RaC')
lam=var('lam') ##lambda
pii=var('pii') ##custom pi variable to avoid memory leaks reported with RDF(pi) calls
##End Declare Variables

##Begin Inputs: Numerical Values
#temporally periodic case
#anum=RDF(100.0); bnum=RDF(100.0)  ##Numerical values of a and b
#lamnum=RDF(1.0); ##aspect ratio
#xnum=RDF(0.01*lamnum); znum=RDF(0.5); tnum=RDF(0.005) ##Numerical value of x,z, and t
#knum=RDF(35.0) ##controls sharpness in logistic function
#zInum=RDF(0.5) ##vertical position of compositional interface
#RaTnum=RDF(1.0e5); RaCnum=RDF(1.0e5) ##Rayleigh numbers
#end temporally periodic case

#approaching steady state case
piinum=3.1415926535897932
anum=RDF(600.0/(piinum*13.0^0.5)); bnum=RDF(100.0)  ##Numerical values of a and b
cnum=RDF(50.0); dnum=RDF(4500.0/(piinum*13.0^0.5))
lamnum=RDF(3.0/2.0); ##aspect ratio
xnum=RDF(0.01*lamnum); znum=RDF(0.5); tnum=RDF(0.05) ##Numerical value of x,z, and t
knum=RDF(35.0) ##controls sharpness in logistic function
zInum=RDF(0.2) ##vertical position of compositional interface
RaTnum=RDF(1.0e6); RaCnum=RDF(8.0e5) ##Rayleigh numbers
#end approaching steady state case
##End Inputs: Numerical Values

####Begin Functions -- use real_part to ensure result is of real type, and to suppress small imaginary parts due to numerical error 
f=function('f')
psi=f(t)*sin(pii*x/lam)*sin(pii*z) ##presumed stream function -- use generic time dependence
#f_trial(t)=a*sin(b*pii*t) ##temporally periodic example
f_trial(t)=a*sin(b*pii*t)*e^(-c*t)+d ##approaching steady state example

D=sin(pii*x/lam)*sin(pii*z)
z_bot=(1/pii)*arcsin(D)
z_top=1-(1/pii)*arcsin(D)

S=2*unit_step(-x+lam/2)-1

#######################################NEED TO UPDATE THIS FORMULA TO USE THE MORE GENERAL TIME DEPENDENCE
tp=var('tp') ##t prime -- variable of integration
#z0=(1/pii)*arccos((jacobi('cn',elliptic_f(pii*z,1/D^2)-S*(I*pii^2*D/lam)*integrate(f(tp),tp,0,t),1/D^2))) ##no real_part
z0=(1/pii)*arccos(real_part(jacobi('cn',elliptic_f(pii*z,1/D^2)-S*(I*pii^2*D/lam)*integrate(f(tp),tp,0,t),1/D^2))) ##real_part used for plotting

##testing
Q=log(abs(csc(pii*z)+cot(pii*z)))-(pii/lam)*S(x=xb)*integrate(f(tp),tp,0,t)
Z0_ad=-1/2*(e^(-Q)-e^Q)
#Z0_p=1/2*(e^(-Q)-e^Q)
#Z0_m=-1/2*(e^(-Q)-e^Q)
#z0_sides_p=(1/pii)*arccot(Z0_p)
#z0_sides_m=(1/pii)*arccot(Z0_m)
z0_sides_Z0_p=(1/pii)*arccot(Z0_ad)
z0_sides_Z0_n=1+(1/pii)*arccot(Z0_ad)
##end testing

C=(1+exp(-2*k*(zI-z0)))^(-1)
####End Functions

###Begin Calculate RMS velocity
u=diff(psi,z) ##velocity components
w=-diff(psi,x)
assume(lam>0)
vrms=sqrt(integrate(integrate(u^2+w^2,x,0,lam),z,0,1)/lam/1) #compute RMS velocity
print("u=",u)
print("LaTeX:", latex(u))
print("w=",w)
print("LaTeX:", latex(w))
print("")
print("vRMS=",vrms)
print("LaTeX:", latex(vrms.factor()))
print("")
print("tnum=",tnum)
print("vRMS(tnum)=",vrms.substitute_function(f,f_trial)(a=anum,b=bnum,c=cnum,d=dnum,t=tnum,lam=lamnum,pii=piinum).n()) ##print numerical value of RMS velocity at t=tnum
print("")
vrms_paper=(pii*sqrt(lam^2+1)/2/lam)*abs(f(t)) ##form in the manuscript after manual simplification
print("vRMS_paper(tnum)=",vrms_paper.substitute_function(f,f_trial)(a=anum,b=bnum,c=cnum,d=dnum,t=tnum,lam=lamnum,pii=piinum).n()) ##print numerical value of RMS velocity at t=tnum
print("")
vrms_pos=(pii*sqrt(lam^2+1)/2/lam)*f_trial(tp) ##this form assumes RMS velocity is positive
print(vrms_pos)
transits_pos=integrate(vrms_pos,tp,0,t,algorithm="giac")
print("transits_pos=",transits_pos)
print("transits_pos(tnum)=",transits_pos(a=anum,b=bnum,c=cnum,d=dnum,t=tnum,lam=lamnum,pii=piinum).n())
print("")
###End Calculate RMS velocity

########Compute Entrainment
#E=integrate(integrate(C,x,0,lam),z,zR,1)/(lam*zI) ## closed form expression may not exist for this double integral
#print("E=",E)
#print("")
########End Compute Entrainment

####Begin Compute G by Integrating the Biharmonic Equation
eta=1 ##assume constant viscosity
inner1=eta*(diff(psi,x,2)-diff(psi,z,2))
LHS1=diff(inner1,x,2)-diff(inner1,z,2)
inner2=eta*diff(diff(psi,z,1),x,1)
LHS2=4*diff(diff(inner2,z,1),x,1)
LHS=LHS1+LHS2
print("Biharmonic LHS=",LHS)
print("")
fG=(RaT-RaC)*(1-z) ##function of integration (constant in x) -- solve for this
G=integrate(LHS,x,algorithm="giac")+fG
print("G=",G.factor())
print("LaTeX:",latex(G.factor()))
print("")
####End Compute G by Integrating the Biharmonic Equation

######Solve for T
T=(G+RaC*C)/RaT
##dTdx=diff(T,x,1) ###hangs
#print("test:",diff(z0,t,1)) ##works (without real_part in z0)
#print("test:",diff(T,t,1)) ##works and is fast
#print("dTdx=",dTdx)

T_centre=(G+RaC*(1+exp(-2*k*(zI-z)))^(-1))/RaT ##T at (x,z)=(lam/2,1/2)
######End Solve for T

##Test Output
print("x,z,t=",xnum,znum,tnum)
#print("D=",D(x=xnum,z=znum,lam=lamnum).n())
#print("z_top=",z_top(x=xnum,z=znum,lam=lamnum).n())
#print("z_bot=",z_bot(x=xnum,z=znum,lam=lamnum).n())
#print("z0=",z0(x=xnum,z=znum,t=tnum,a=anum,b=bnum,lam=lamnum).n())
#print("z0_sides_p=",z0_sides_p(z=znum,t=tnum,a=anum,b=bnum,lam=lamnum).n())
#print("z0_sides_m=",z0_sides_m(z=znum,t=tnum,a=anum,b=bnum,lam=lamnum).n())
#print("C=",C(x=xnum,z=znum,t=tnum,a=anum,b=bnum,lam=lamnum,k=knum,zI=zInum).n())
#print("T=",T(x=xnum,z=znum,t=tnum,a=anum,b=bnum,lam=lamnum,k=knum,zI=zInum,RaT=RaTnum,RaC=RaCnum).n())

print("z0=",z0.substitute_function(f,f_trial)(x=xnum,z=znum,t=tnum,a=anum,b=bnum,c=cnum,d=dnum,lam=lamnum,pii=piinum).n())

xbnum=RDF(0.0) ##left sidewall
print("xb=",xbnum)
print("Z0_ad=",Z0_ad.substitute_function(f,f_trial)(xb=xbnum,z=znum,t=tnum,a=anum,b=bnum,c=cnum,d=dnum,lam=lamnum,pii=piinum).n())
print("z0_sides_Z0_p=",z0_sides_Z0_p.substitute_function(f,f_trial)(xb=xbnum,z=znum,t=tnum,a=anum,b=bnum,c=cnum,d=dnum,lam=lamnum,pii=piinum).n())
print("z0_sides_Z0_n=",z0_sides_Z0_n.substitute_function(f,f_trial)(xb=xbnum,z=znum,t=tnum,a=anum,b=bnum,c=cnum,d=dnum,lam=lamnum,pii=piinum).n())

xbnum=lamnum ##right sidewall
print("xb=",xbnum)
print("Z0_ad=",Z0_ad.substitute_function(f,f_trial)(xb=xbnum,z=znum,t=tnum,a=anum,b=bnum,c=cnum,d=dnum,lam=lamnum,pii=piinum).n())
print("z0_sides_Z0_p=",z0_sides_Z0_p.substitute_function(f,f_trial)(xb=xbnum,z=znum,t=tnum,a=anum,b=bnum,c=cnum,d=dnum,lam=lamnum,pii=piinum).n())
print("z0_sides_Z0_n=",z0_sides_Z0_n.substitute_function(f,f_trial)(xb=xbnum,z=znum,t=tnum,a=anum,b=bnum,c=cnum,d=dnum,lam=lamnum,pii=piinum).n())
#print("")
##End Test Output

##fill plot array
def fill_plot_array(Nxf,Nzf,x1f,x2f,z1f,z2f,plot_functionf): ##end variables with "f" for "fill"
        #plot_function_fast=fast_callable(plot_functionf,vars=[x,z],domain=CDF)
        plot_function_fast=plot_functionf
        plot_array=np.zeros((Nzf+1,Nxf+1)) ##initialize array with zeros
        for ii in range(Nzf+1): ##first array index corresponds to z values
                zval=float(ii)*((z2f-z1f)/float(Nzf))+z1f
                for kk in range(Nxf+1): ##second array index corresponds to x values
                        xval=float(kk)*((x2f-x1f)/float(Nxf))+x1f
                        plot_array[ii,kk]=plot_function_fast(x=xval,z=zval).n()
                        #print("(x,z)=",xval,zval,plot_array[ii,kk])
        return plot_array
##end plot array

####
def display_plot(Nx,Nz,x1,x2,z1,z2,plot_function,title):
        #plot_array=fill_plot_array(Nx,Nz,x1,x2,z1,z2,plot_function)  ##are the arguments here the problem?
        plot_array=fill_plot_array(Nxf=Nx,Nzf=Nz,x1f=x1,x2f=x2,z1f=z1,z2f=z2,plot_functionf=plot_function) ##use "f" to distinguish variables here
        print("min,max=",plot_array.min(),plot_array.max())
        image=plt.imshow(plot_array,interpolation="kaiser",cmap="copper",vmin=0.0,vmax=1.0) ##image to be plotted
        axes=plt.gca() #axes
        axes.set_ylim(axes.get_ylim()[::-1]) ##invert y axis
        #axes.set_xlabel('x',fontsize=16)
        #axes.set_ylabel('z',fontsize=16)
        axes.axes.xaxis.set_ticks([]) ##hide axes ticks
        axes.axes.yaxis.set_ticks([]) ##hide axes ticks
        #plt.title(title,fontsize=20) #hide title for now
        cbar=plt.colorbar(image)
        cbar.ax.tick_params(labelsize=14)
        plt.show()
####

###create plot datafile
def create_plot_datafile(Nx,Nz,x1,x2,z1,z2,plot_function,fname):
        plot_function_fast=plot_function
        with open(fname, 'w') as f:
                for kk in range(Nz+1):
                        zval=float(kk)*((z2-z1)/float(Nz))+z1
                        for ii in range(Nx+1):
                                xval=float(ii)*((x2-x1)/float(Nx))+x1
                                val=plot_function_fast(x=xval,z=zval).n()
                                print(xval,zval,val, file=f)
###Plots
#display_plot(Nx=11,Nz=11,x1=0.01,x2=0.99,z1=0.01,z2=0.99,plot_function=C(t=0.01,a=anum,b=bnum,k=35.0,zI=0.5,lam=1.0),title="C")
#display_plot(Nx=121,Nz=121,x1=0.001,x2=0.999,z1=0.001,z2=0.999,plot_function=T(t=0.01,a=anum,b=bnum,k=35.0,zI=0.5,lam=1.0,RaT=1.e5,RaC=1.e5),title="T") ##for figures
#display_plot(Nx=181,Nz=121,x1=0.001,x2=1.4999,z1=0.001,z2=0.999,plot_function=C.substitute_function(f,f_trial)(t=0.05,a=anum,b=bnum,c=cnum,d=dnum,k=knum,zI=zInum,lam=lamnum,RaT=RaTnum,RaC=RaCnum),title="C") ##steady state case

##example for steady state T datafile
#T_in=T.substitute_function(f,f_trial)(t=0.05,a=anum,b=bnum,c=cnum,d=dnum,k=knum,zI=zInum,lam=lamnum,pii=piinum,RaT=RaTnum,RaC=RaCnum)
#create_plot_datafile(Nx=451,Nz=301,x1=0.001,x2=1.4999,z1=0.001,z2=0.999,plot_function=T_in,fname='T_steady_state_451x301_t0p05.dat')

###End Plots
