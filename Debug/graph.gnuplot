set autoscale
set xtic auto
set ytic auto

set title "|psi|^2"
set xlabel "x [Hartree]"
set ylabel "|psi|^2"
#set xr [-3.15:3.15]
#set yr [-1.5:1.5]

set style fill transparent solid 0.75 noborder

plot  "output.dat" using 1:4 with lines lt 1 title "|psi|^2",\
      ""           using 1:2 with filledcurves x1 lt 2 title "V(x)"