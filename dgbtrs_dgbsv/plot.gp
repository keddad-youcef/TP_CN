set title "Resolution Methods Performance"
set xlabel "Taille de Matrice"
set ylabel "Temps d'execution (sec)"
set auto x
plot 'Performance.dat' using 1:2 title 'DGBTRS' with lines, 'Performance.dat' using 1:3 title 'DGBSV' with lines
pause -1
