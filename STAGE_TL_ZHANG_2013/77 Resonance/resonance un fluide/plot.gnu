plot "ETA_Goutte_Attachee_Inviscid.dat" w l,  "ETA_Goutte_Attachee_Inviscid.dat" u 1:3 w l 
pause 15
set term postscript
set output "ETA_Goutte_Attachee_Inviscid.eps" 
replot
quit
