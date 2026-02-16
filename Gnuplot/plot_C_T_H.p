#gnuplot: plot C,T, and H snapshots
reset

#### Input Parameters

lam=1.0 ##aspect ratio

## Size of plot windows
term_x=640 ##terminal size along x -- default 640 (850 for problem 2)
term_y=480 ##terminal size along y -- default 480

## Size of data points on plots
psize=0.125 ##default is 0.125 but may be decreased if data points overlap or increased if white space appears between data points

##file names for C, T and H
fnameC="../Fortran/C_data.dat"
fnameT="../Fortran/T_data.dat"
fnameH="../Fortran/H_data.dat"

# Internal Variables
x1=0.0; x2=lam
z1=0.0; z2=1.0

##get min and max values
stats fnameC u 3; Cmin=STATS_min; Cmax=STATS_max
stats fnameT u 3; Tmin=STATS_min; Tmax=STATS_max
stats fnameH u 3; Hmin=STATS_min; Hmax=STATS_max

##options for all plots
set lmargin 0; set rmargin 0; set tmargin 0; set bmargin 0
set xrange [x1:x2]; set yrange [z1:z2];
unset xtics; unset ytics
set cbtics out; set cbtics nomirror; set cbtics offset -1
set size ratio -1
set view map
unset key
unset border
load 'batlow.pal'

##C plot
set term qt 0 enhanced font "Times,20" title "C" size term_x,term_y
set cbrange [Cmin:Cmax]
set cbtics ( Cmin, 0.25, 0.5, 0.75, Cmax )
set format cb " %.2f"
splot fnameC u 1:2:3 w p pt 5 ps psize palette

##T plot
set term qt 1 enhanced font "Times,20" title "T" size term_x,term_y
set cbrange [Tmin:Tmax]
set cbtics ( 0, 0.25, 0.5, 0.75, 1 ) 
#set cbtics ( Tmin, 0.25, 0.5, 0.75, Tmax ) #t=0
#set cbtics ( Tmin, 0, 0.25, 0.5, 0.75, 1, Tmax ) #general
set format cb " %.2f"
splot fnameT u 1:2:3 w p pt 5 ps psize palette

##H plot
set term qt 2 enhanced font "Times,20" title "H" size term_x,term_y
######symmetric log scale and inverse
symlog(z)  = (-1.0 < z && z < 1.0) ? z/10. \
           : (z < 0.0) ? -log10(-z) - 0.1 \
           : log10(z) + 0.1

invsymlog(z) = (-0.1 < z && z < 0.1) ? z*10. \
             : (z < 0.0) ? -10.0**(-(z+0.1)) \
             : 10.0**(z-0.1)

set nonlinear cb via symlog(z) inv invsymlog(z) ##axis is "cb" and dummy variable is "z" -- see gnuplot help
######end symmetric log scale and inverse
set cbrange [Hmin:Hmax]
set cbtics ( Hmin, -1e10, -1e8, -1e6,-1e4,-1e2, 0, 1e2, 1e4, 1e6, 1e8, 1e10, Hmax ) # tics out-of-range will not appear
set format cb " % .1t x 10^{%T}"
splot fnameH u 1:2:3 w p pt 5 ps psize palette
