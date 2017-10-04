plot "ETA_SPHERE.dat" w l,  "ETA_SPHERE.dat" u 1:3 w l 
pause 5
set term pdf
set output "ETA_SPHERE.pdf" 
replot
quit
