set autoscale
set xtic auto
set ytic auto

set grid

set title "Infinite Square Well"
set xlabel "x [Hartree]"
set ylabel "|psi|^2"
#set xr [-3.15:3.15]
#set yr [-1.5:1.5]

set style fill transparent solid 0.75 noborder

plot  "isw_n1.dat" using 1:4 with lines lt 1 title "n = 0",\
      "isw_n2.dat" using 1:4 with lines lt 2 title "n = 1",\
      "isw_n3.dat" using 1:4 with lines lt 3 title "n = 2",\
      "isw_n4.dat" using 1:4 with lines lt 4 title "n = 3",\
      "isw_n5.dat" using 1:4 with lines lt 5 title "n = 4",\
      "isw_n6.dat" using 1:4 with lines lt 6 title "n = 5",\
      "isw_n7.dat" using 1:4 with lines lt 7 title "n = 6",\
      ""           using 1:2 with filledcurves x1 lt 8 title "V(x)"