      # Gnuplot script file for plotting data in file "cmake-build-debug-wsl/output.dat"
      # Plot with:
      # terminal```
      # >gnuplot
      # >load 'graphs.gnuplot'
      # ```
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set title "Ising Model"
      set xlabel "Temperature []"
      set ylabel "Magnetization []"
#      set xr [-3.15:3.15]
#      set yr [-1.5:1.5]
      plot    "cmake-build-debug-wsl/output.dat" using 1:4 with points