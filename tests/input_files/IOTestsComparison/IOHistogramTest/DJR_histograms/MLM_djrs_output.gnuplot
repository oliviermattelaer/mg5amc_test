
################################################################################
#
# This gnuplot file was generated by MadGraph5_aMC@NLO project, a program which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond. It also perform the
# integration and/or generate events for these processes, at LO and NLO accuracy.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################
# Automatic plotting from MG5aMC
reset

set lmargin 10
set rmargin 0
set terminal postscript portrait enhanced mono dashed lw 1.0 "Helvetica" 9 
# The pdf terminal offers transparency support, but you will have to adapt things a bit
#set terminal pdf enhanced font "Helvetica 12" lw 1.0 dashed size 29.7cm, 21cm
set key font ",9"
set key samplen "2"
set output "MLM_djrs_output.ps"

# This is the "PODO" color palette of gnuplot v.5, but with the order
# changed: palette of colors selected to be easily distinguishable by
# color-blind individuals with either protanopia or deuteranopia. Bang
# Wong [2011] Nature Methods 8, 441.

set style line  1 lt 1 lc rgb "#009e73" lw 2.5
set style line 11 lt 2 lc rgb "#009e73" lw 2.5
set style line 21 lt 4 lc rgb "#009e73" lw 2.5
set style line 31 lt 6 lc rgb "#009e73" lw 2.5
set style line 41 lt 8 lc rgb "#009e73" lw 2.5

set style line  2 lt 1 lc rgb "#0072b2" lw 2.5
set style line 12 lt 2 lc rgb "#0072b2" lw 2.5
set style line 22 lt 4 lc rgb "#0072b2" lw 2.5
set style line 32 lt 6 lc rgb "#0072b2" lw 2.5
set style line 42 lt 8 lc rgb "#0072b2" lw 2.5

set style line  3 lt 1 lc rgb "#d55e00" lw 2.5
set style line 13 lt 2 lc rgb "#d55e00" lw 2.5
set style line 23 lt 4 lc rgb "#d55e00" lw 2.5
set style line 33 lt 6 lc rgb "#d55e00" lw 2.5
set style line 43 lt 8 lc rgb "#d55e00" lw 2.5

set style line  4 lt 1 lc rgb "#f0e442" lw 2.5
set style line 14 lt 2 lc rgb "#f0e442" lw 2.5
set style line 24 lt 4 lc rgb "#f0e442" lw 2.5
set style line 34 lt 6 lc rgb "#f0e442" lw 2.5
set style line 44 lt 8 lc rgb "#f0e442" lw 2.5

set style line  5 lt 1 lc rgb "#56b4e9" lw 2.5
set style line 15 lt 2 lc rgb "#56b4e9" lw 2.5
set style line 25 lt 4 lc rgb "#56b4e9" lw 2.5
set style line 35 lt 6 lc rgb "#56b4e9" lw 2.5
set style line 45 lt 8 lc rgb "#56b4e9" lw 2.5

set style line  6 lt 1 lc rgb "#cc79a7" lw 2.5
set style line 16 lt 2 lc rgb "#cc79a7" lw 2.5
set style line 26 lt 4 lc rgb "#cc79a7" lw 2.5
set style line 36 lt 6 lc rgb "#cc79a7" lw 2.5
set style line 46 lt 8 lc rgb "#cc79a7" lw 2.5

set style line  7 lt 1 lc rgb "#e69f00" lw 2.5
set style line 17 lt 2 lc rgb "#e69f00" lw 2.5
set style line 27 lt 4 lc rgb "#e69f00" lw 2.5
set style line 37 lt 6 lc rgb "#e69f00" lw 2.5
set style line 47 lt 8 lc rgb "#e69f00" lw 2.5

set style line  8 lt 1 lc rgb "black" lw 2.5
set style line 18 lt 2 lc rgb "black" lw 2.5
set style line 28 lt 4 lc rgb "black" lw 2.5
set style line 38 lt 6 lc rgb "black" lw 2.5
set style line 48 lt 7 lc rgb "black" lw 2.5


set style line 999 lt 1 lc rgb "gray" lw 2.5

safe(x,y,a) = (y == 0.0 ? a : x/y)

set style data histeps



################################################################################
### Rendering of the plot titled 'log10d01'
################################################################################

set multiplot
set label "log10d01" font ",13" at graph 0.04, graph 1.05
set xrange [0.0000e+00:2.8200e+00]
set bmargin 0 
set tmargin 0
set xtics nomirror
set ytics nomirror
set mytics 10
set xtics auto
set key horizontal noreverse maxcols 1 width -4 
set label front 'MadGraph5\_aMC\@NLO' font "Courier,11" rotate by 90 at graph 1.02, graph 0.04

#-- rendering subhistograms 'None and None results'

set format y '10^{%%T}'
set yrange [1.3539e-15:3.7328e-11]
set origin 0.0000e+00, 5.0000e-01
set size 1.0000e+00, 4.0000e-01
set mytics 10
set ytics auto
set format x ''
set logscale y
set ylabel "{/Symbol s} per bin [pb]"

plot \
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):13 ls 31 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):12 ls 31 title 'all jet samples, merging scale variation',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):7 ls 11 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):6 ls 11 title 'all jet samples, scale variation',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):10 ls 21 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):9 ls 21 title 'all jet samples, PDF variation',\
'MLM_djrs_output.HwU' index 3 using (($1+$2)/2):3:4 w yerrorbar ls 4 title '',\
'MLM_djrs_output.HwU' index 3 using (($1+$2)/2):3 ls 4 title 'jet sample 2',\
'MLM_djrs_output.HwU' index 2 using (($1+$2)/2):3:4 w yerrorbar ls 3 title '',\
'MLM_djrs_output.HwU' index 2 using (($1+$2)/2):3 ls 3 title 'jet sample 1',\
'MLM_djrs_output.HwU' index 1 using (($1+$2)/2):3:4 w yerrorbar ls 2 title '',\
'MLM_djrs_output.HwU' index 1 using (($1+$2)/2):3 ls 2 title 'jet sample 0',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):3:4 w yerrorbar ls 1 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):3 ls 1 title 'log10d01, all jet samples'
#-- rendering subhistograms 'Relative scale and PDF uncertainty'
unset label
unset format
set yrange [-1.9800e-01:2.1396e-01]
set origin 0.0000e+00, 3.5000e-01
set size 1.0000e+00, 1.5000e-01
set mytics 2
set ytics auto
set format x
unset logscale y
set ylabel "(1) rel.unc."
set label "Relative uncertainties w.r.t. central value(s)" font ",9" front at graph 0.03, graph 0.13
plot \
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):(safe($13,$3,1.0)-1.0) ls 31 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):(safe($12,$3,1.0)-1.0) ls 31 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):(safe($7,$3,1.0)-1.0) ls 11 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):(safe($6,$3,1.0)-1.0) ls 11 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):(safe($10,$3,1.0)-1.0) ls 21 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):(safe($9,$3,1.0)-1.0) ls 21 title '',\
0.0 ls 999 title '',\
'MLM_djrs_output.HwU' index 0 using (($1+$2)/2):(0.0):(safe($4,$3,0.0)) w yerrorbar ls 1 title ''

unset label

################################################################################

################################################################################
### Rendering of the plot titled 'log10d12'
################################################################################

set multiplot
set label "log10d12" font ",13" at graph 0.04, graph 1.05
set xrange [0.0000e+00:2.6100e+00]
set bmargin 0 
set tmargin 0
set xtics nomirror
set ytics nomirror
set mytics 10
set xtics auto
set key horizontal noreverse maxcols 1 width -4 
set label front 'MadGraph5\_aMC\@NLO' font "Courier,11" rotate by 90 at graph 1.02, graph 0.04

#-- rendering subhistograms 'None and None results'

set format y '10^{%%T}'
set yrange [1.1712e-15:2.7041e-11]
set origin 0.0000e+00, 5.0000e-01
set size 1.0000e+00, 4.0000e-01
set mytics 10
set ytics auto
set format x ''
set logscale y
set ylabel "{/Symbol s} per bin [pb]"

plot \
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):13 ls 31 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):12 ls 31 title 'all jet samples, merging scale variation',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):7 ls 11 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):6 ls 11 title 'all jet samples, scale variation',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):10 ls 21 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):9 ls 21 title 'all jet samples, PDF variation',\
'MLM_djrs_output.HwU' index 7 using (($1+$2)/2):3:4 w yerrorbar ls 4 title '',\
'MLM_djrs_output.HwU' index 7 using (($1+$2)/2):3 ls 4 title 'jet sample 2',\
'MLM_djrs_output.HwU' index 6 using (($1+$2)/2):3:4 w yerrorbar ls 3 title '',\
'MLM_djrs_output.HwU' index 6 using (($1+$2)/2):3 ls 3 title 'jet sample 1',\
'MLM_djrs_output.HwU' index 5 using (($1+$2)/2):3:4 w yerrorbar ls 2 title '',\
'MLM_djrs_output.HwU' index 5 using (($1+$2)/2):3 ls 2 title 'jet sample 0',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):3:4 w yerrorbar ls 1 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):3 ls 1 title 'log10d12, all jet samples'
#-- rendering subhistograms 'Relative scale and PDF uncertainty'
unset label
unset format
set yrange [-1.9551e-01:2.0762e-01]
set origin 0.0000e+00, 3.5000e-01
set size 1.0000e+00, 1.5000e-01
set mytics 2
set ytics auto
set format x
unset logscale y
set ylabel "(1) rel.unc."
set label "Relative uncertainties w.r.t. central value(s)" font ",9" front at graph 0.03, graph 0.13
plot \
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):(safe($13,$3,1.0)-1.0) ls 31 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):(safe($12,$3,1.0)-1.0) ls 31 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):(safe($7,$3,1.0)-1.0) ls 11 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):(safe($6,$3,1.0)-1.0) ls 11 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):(safe($10,$3,1.0)-1.0) ls 21 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):(safe($9,$3,1.0)-1.0) ls 21 title '',\
0.0 ls 999 title '',\
'MLM_djrs_output.HwU' index 4 using (($1+$2)/2):(0.0):(safe($4,$3,0.0)) w yerrorbar ls 1 title ''

unset label

################################################################################

################################################################################
### Rendering of the plot titled 'log10d23'
################################################################################

set multiplot
set label "log10d23" font ",13" at graph 0.04, graph 1.05
set xrange [0.0000e+00:2.4300e+00]
set bmargin 0 
set tmargin 0
set xtics nomirror
set ytics nomirror
set mytics 10
set xtics auto
set key horizontal noreverse maxcols 1 width -4 
set label front 'MadGraph5\_aMC\@NLO' font "Courier,11" rotate by 90 at graph 1.02, graph 0.04

#-- rendering subhistograms 'None and None results'

set format y '10^{%%T}'
set yrange [1.1733e-15:2.4993e-11]
set origin 0.0000e+00, 5.0000e-01
set size 1.0000e+00, 4.0000e-01
set mytics 10
set ytics auto
set format x ''
set logscale y
set ylabel "{/Symbol s} per bin [pb]"

plot \
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):13 ls 31 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):12 ls 31 title 'all jet samples, merging scale variation',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):7 ls 11 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):6 ls 11 title 'all jet samples, scale variation',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):10 ls 21 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):9 ls 21 title 'all jet samples, PDF variation',\
'MLM_djrs_output.HwU' index 11 using (($1+$2)/2):3:4 w yerrorbar ls 4 title '',\
'MLM_djrs_output.HwU' index 11 using (($1+$2)/2):3 ls 4 title 'jet sample 2',\
'MLM_djrs_output.HwU' index 10 using (($1+$2)/2):3:4 w yerrorbar ls 3 title '',\
'MLM_djrs_output.HwU' index 10 using (($1+$2)/2):3 ls 3 title 'jet sample 1',\
'MLM_djrs_output.HwU' index 9 using (($1+$2)/2):3:4 w yerrorbar ls 2 title '',\
'MLM_djrs_output.HwU' index 9 using (($1+$2)/2):3 ls 2 title 'jet sample 0',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):3:4 w yerrorbar ls 1 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):3 ls 1 title 'log10d23, all jet samples'
#-- rendering subhistograms 'Relative scale and PDF uncertainty'
unset label
unset format
set yrange [-2.1261e-01:2.3059e-01]
set origin 0.0000e+00, 3.5000e-01
set size 1.0000e+00, 1.5000e-01
set mytics 2
set ytics auto
set format x
unset logscale y
set ylabel "(1) rel.unc."
set label "Relative uncertainties w.r.t. central value(s)" font ",9" front at graph 0.03, graph 0.13
plot \
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):(safe($13,$3,1.0)-1.0) ls 31 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):(safe($12,$3,1.0)-1.0) ls 31 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):(safe($7,$3,1.0)-1.0) ls 11 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):(safe($6,$3,1.0)-1.0) ls 11 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):(safe($10,$3,1.0)-1.0) ls 21 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):(safe($9,$3,1.0)-1.0) ls 21 title '',\
0.0 ls 999 title '',\
'MLM_djrs_output.HwU' index 8 using (($1+$2)/2):(0.0):(safe($4,$3,0.0)) w yerrorbar ls 1 title ''

unset label

################################################################################

################################################################################
### Rendering of the plot titled 'log10d34'
################################################################################

set multiplot
set label "log10d34" font ",13" at graph 0.04, graph 1.05
set xrange [0.0000e+00:2.3400e+00]
set bmargin 0 
set tmargin 0
set xtics nomirror
set ytics nomirror
set mytics 10
set xtics auto
set key horizontal noreverse maxcols 1 width -4 
set label front 'MadGraph5\_aMC\@NLO' font "Courier,11" rotate by 90 at graph 1.02, graph 0.04

#-- rendering subhistograms 'None and None results'

set format y '10^{%%T}'
set yrange [1.1525e-15:2.6431e-11]
set origin 0.0000e+00, 5.0000e-01
set size 1.0000e+00, 4.0000e-01
set mytics 10
set ytics auto
set format x ''
set logscale y
set ylabel "{/Symbol s} per bin [pb]"

plot \
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):13 ls 31 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):12 ls 31 title 'all jet samples, merging scale variation',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):7 ls 11 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):6 ls 11 title 'all jet samples, scale variation',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):10 ls 21 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):9 ls 21 title 'all jet samples, PDF variation',\
'MLM_djrs_output.HwU' index 15 using (($1+$2)/2):3:4 w yerrorbar ls 4 title '',\
'MLM_djrs_output.HwU' index 15 using (($1+$2)/2):3 ls 4 title 'jet sample 2',\
'MLM_djrs_output.HwU' index 14 using (($1+$2)/2):3:4 w yerrorbar ls 3 title '',\
'MLM_djrs_output.HwU' index 14 using (($1+$2)/2):3 ls 3 title 'jet sample 1',\
'MLM_djrs_output.HwU' index 13 using (($1+$2)/2):3:4 w yerrorbar ls 2 title '',\
'MLM_djrs_output.HwU' index 13 using (($1+$2)/2):3 ls 2 title 'jet sample 0',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):3:4 w yerrorbar ls 1 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):3 ls 1 title 'log10d34, all jet samples'
#-- rendering subhistograms 'Relative scale and PDF uncertainty'
unset label
unset format
set yrange [-2.6022e-01:2.8412e-01]
set origin 0.0000e+00, 3.5000e-01
set size 1.0000e+00, 1.5000e-01
set mytics 2
set ytics auto
set format x
unset logscale y
set ylabel "(1) rel.unc."
set label "Relative uncertainties w.r.t. central value(s)" font ",9" front at graph 0.03, graph 0.13
plot \
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):(safe($13,$3,1.0)-1.0) ls 31 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):(safe($12,$3,1.0)-1.0) ls 31 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):(safe($7,$3,1.0)-1.0) ls 11 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):(safe($6,$3,1.0)-1.0) ls 11 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):(safe($10,$3,1.0)-1.0) ls 21 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):(safe($9,$3,1.0)-1.0) ls 21 title '',\
0.0 ls 999 title '',\
'MLM_djrs_output.HwU' index 12 using (($1+$2)/2):(0.0):(safe($4,$3,0.0)) w yerrorbar ls 1 title ''

unset label

################################################################################
unset multiplot
!ps2pdf "MLM_djrs_output.ps" &> /dev/null