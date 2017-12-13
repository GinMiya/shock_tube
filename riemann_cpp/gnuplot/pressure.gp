reset
set xr [0:1]
set yr [0:1.1]
set ticslevel 0
set size square

set grid
set term gif animate delay 20
set output "pressure.gif"
unset title
set tics font "Arial,20"
set xl 'x' font "Arial,20"
set yl 'pressure' font "Arial,20"

n0 = 0
nm = 20
dn = 1

load './gnuplot/pressure_gif.plt'
