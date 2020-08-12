set autoscale
set xtic auto
set ytic auto

set grid

set title "Quantum Harmonic Oscillator"
set xlabel "x [Hartree]"
set ylabel "|psi|^2"

set style fill transparent solid 0.75 noborder

plot  "qho_n0.dat" using 1:4 with lines lt 1 title "n = 0",\
      "qho_n1.dat" using 1:4 with lines lt 2 title "n = 1",\
      "qho_n2.dat" using 1:4 with lines lt 3 title "n = 2",\
      "qho_n3.dat" using 1:4 with lines lt 4 title "n = 3",\
      "qho_n4.dat" using 1:4 with lines lt 5 title "n = 4",\
      "qho_n5.dat" using 1:4 with lines lt 6 title "n = 5",\
      "qho_n6.dat" using 1:4 with lines lt 7 title "n = 6",\
      ""           using 1:2 with filledcurves x1 lt 8 title "V(x)"