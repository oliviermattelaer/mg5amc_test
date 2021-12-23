
# General settings:
   set output '<name>.eps'
   set terminal postscript eps enhanced color 'Helvetica' 16
   set style data lines
   set size 1.0,2.0

   set multiplot

   reset
   set xlabel "p_{{/Symbol \\136} evol}  (GeV)"

   #set ylabel "Deviation [%]"
   #set format y "    %1.1f"
   set logscale y

   set key right top

   set origin 0.,0.

plot "<file1>" u 1:2 w l, "<file2>" u 1:2 w l


  set origin 0.,1.
plot "<file3>" u 1:2 w l, "<file4>" u 1:2 w l
