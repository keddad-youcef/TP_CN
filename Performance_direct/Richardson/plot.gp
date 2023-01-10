set title "Methode de Richarson"
set xlabel "Alpha"
set ylabel "Nombre d'iterations pour la convergence"
set grid
set style data linespoints
plot 'Alpha_Richarson.dat' using 1:2 lt 1 lw 2 smooth bezier title 'Alpha' with lines
pause -1
