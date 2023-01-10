reset
set title "Resolution Methods Performance"
set xlabel "Taille de Matrice"
set ylabel "Temps d'execution (sec)"
set auto x
set grid
set style data linespoints
plot 'Performance LU.dat' using 1:2 lt 1 lw 2 smooth bezier title 'LU_FACTO', 'Performance LU.dat' using 1:3 lt 2 lw 2 smooth bezier title 'DGBTRF'
pause -1
