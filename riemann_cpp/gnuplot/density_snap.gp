reset
set xr [0:1]
set yr [0:1.1]
set ticslevel 0
set size square

set grid
set term png
set output "./data/density_snap.png"
set tics font "Arial,20"
set xl 'x' font "Arial,20"
set yl 'density' font "Arial,20"

pl './data/shock_tube_19.dat' u 1:2 w p pt 7 ps 0.5 notitle
