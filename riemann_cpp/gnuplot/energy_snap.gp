reset
set xr [0:1]
set yr [1.5:3]
set ticslevel 0
set size square

set grid
set term png
set output "./data/energy_snap.png"
set tics font "Arial,20"
set xl 'x' font "Arial,20"
set yl 'energy' font "Arial,20"

pl './data/shock_tube_19.dat' u 1:5 w p pt 7 ps 0.5 notitle
