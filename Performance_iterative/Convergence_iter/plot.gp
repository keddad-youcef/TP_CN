set title "Convergence history iterative"
set xlabel "Nombre d'teration"
set ylabel "l'erreur residuel"
set auto x
set grid
set style data linespoints
plot 'Gauss_Seidel.dat' using 1:2 lt 3 lw 3 smooth bezier title 'Gauss Seidel'
pause -1
