set title "Convergence history iterative"
set xlabel "Nombre d'teration"
set ylabel "l'erreur residuel"
set auto x
set grid
set style data linespoints
plot 'Richardson.dat' using 1:2 lt 1 lw 2 smooth bezier title 'Richardson', 'Jacobi.dat' using 1:2 lt 2 lw 2 smooth bezier title 'Jacobi', 'Gauss_Seidel.dat' using 1:2 lt 3 lw 3 smooth bezier title 'Gauss Seidel'
pause -1
