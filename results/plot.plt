# Plot grid-precision data
filename = "/home/boris/Documenti/grid-precision/results/spline00.dat";

set terminal wxt 0 title 'x'
set xlabel "Time"
set ylabel "x"
set grid
plot filename using 1:2 w l title 'x'




