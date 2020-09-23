set terminal pdfcairo
set output "Punto2B.pdf"
set xrange [0.1:15.0]
set yrange [-0.5:1.0]
set xlabel "l"
set ylabel "R(r = 1, l)"
set title "R(r=1,l) - Function For Finding Eigenvalues"
plot "Punto2B.txt" u 1:2 w l lw 2.5 t "Computed Data"
