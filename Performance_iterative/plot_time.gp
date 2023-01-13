set title "Time Perf for iterative methods"
set xlabel "taille de la matrice"
set ylabel "Temps d'execution (s)"
set auto x
set grid
set style data linespoints
plot 'perf_time.dat' using 1:2 lt 1 lw 2 smooth bezier title 'Jacobi', 'perf_time.dat' using 1:3 lt 2 lw 2 smooth bezier title 'Gauss Seidel'
pause -1
