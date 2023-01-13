set title "Convergence en fonction de la taille de matrice"
set xlabel "Tolerence (10^-)"
set ylabel "nombre d'iterations"
set auto x
set grid
set style data linespoints
plot 'perf_tol.dat' using 1:2 lt 1 lw 2 smooth bezier title 'Richardson', 'perf_tol.dat' using 1:3 lt 2 lw 2 smooth bezier title 'Jacobi', 'perf_tol.dat' using 1:4 lt 3 lw 3 smooth bezier title 'Gauss Seidel'
pause -1
