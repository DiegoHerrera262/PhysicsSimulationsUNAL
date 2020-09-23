set terminal pdfcairo
set output "Punto2A.pdf"
set xrange [0.0:10.0]
set yrange [-0.5:1.0]
set xlabel "r"
set ylabel "R(r)"
set title "Radial Eigen Function for Drum Normal Mode lambda = 1.0"
plot "Punto2A.txt" u 2:3 w l lw 2.5 t "Lambda = 1.0"
