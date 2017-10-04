plot "MODE_SPHERE.dat" w l,  "MODE_SPHERE.dat" u 1:4 w l 
pause 5
set term pdf
set output "MODE_SPHERE.pdf" 
replot
quit
