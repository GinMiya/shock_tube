if(exist ("n") ==0 || n < 0 ) n = n0

filename = sprintf("./data/shock_tube_%d.dat",n)
set title sprintf("t = %f",n*0.025)
plot filename using 1:5 w l notitle

n = n + dn
if ( n < nm ) reread
undefine n
