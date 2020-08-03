set autoscale
set xtic auto
set ytic auto

set title "QHO"
set xlabel "x [Hartree]"
set ylabel "|psi|^2"
#set xr [-3.15:3.15]
#set yr [-1.5:1.5]

set style fill transparent solid 0.5 noborder

plot  "n0.dat" using 1:2 with filledcurves x1 lt 1 title "V(x)",\
      "n0.dat" using 1:4 with lines lw 3 lt 2 title "n = 0",\
      "n1.dat" using 1:4 with lines lw 3 lt 3 title "n = 1",\
      "n2.dat" using 1:4 with lines lw 3 lt 4 title "n = 2",\
      "n3.dat" using 1:4 with lines lw 3 lt 5 title "n = 3",\
      "n4.dat" using 1:4 with lines lw 3 lt 6 title "n = 4",\
      "n5.dat" using 1:4 with lines lw 3 lt 7 title "n = 5",\
      "n6.dat" using 1:4 with lines lw 3 lt 8 title "n = 6",\