#format output file for PDALAC postprocessing
m="./data.dat"
set nokey
set grid
set title 'Stress vs Strain'
set term wxt 0
set xlabel "Strain(1)"
set ylabel "Stress(1)"
set autoscale xy
plot m using 1:7 with linespoints
set autoscale xy
set term wxt 1
set xlabel "Strain(2)"
set ylabel "Stress(2)"
set autoscale xy
plot m using 2:8 with linespoints
set autoscale xy
set term wxt 2
set xlabel "Strain(3)"
set ylabel "Stress(3)"
set autoscale xy
plot m using 3:9 with linespoints
set autoscale xy
set term wxt 3
