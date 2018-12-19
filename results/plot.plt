# Plot grid-precision data
# filename = "/home/boris/Documenti/grid-precision/results/spline00.dat";

# set terminal wxt 0 title 'x'
# set xlabel "Time"
# set ylabel "x"
# set grid
# plot filename using 1:2 w l title 'x'

filename = "/home/boris/Documenti/grid-precision/results/grid_precision.dat";

set terminal wxt 0 title 'x'
set style data histogram
set xlabel "Time"
set ylabel "x"
set grid
set logscale y
plot filename using 2 title 'x'

set terminal wxt 1 title 'y'
set style data histogram
set xlabel "Time"
set ylabel "y"
set grid
set logscale y
plot filename using 3 title 'y'

set terminal wxt 2 title 'z'
set style data histogram
set xlabel "Time"
set ylabel "z"
set grid
set logscale y
plot filename using 4 title 'z'

set terminal wxt 3 title 'vx'
set style data histogram
set xlabel "Time"
set ylabel "vx"
set grid
set logscale y
plot filename using 5 title 'vx'

set terminal wxt 4 title 'vy'
set style data histogram
set xlabel "Time"
set ylabel "vy"
set grid
set logscale y
plot filename using 6 title 'vy'

set terminal wxt 5 title 'vz'
set style data histogram
set xlabel "Time"
set ylabel "vz"
set grid
set logscale y
plot filename using 7 title 'vz'







