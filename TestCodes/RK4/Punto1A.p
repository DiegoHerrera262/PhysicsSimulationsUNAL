set terminal pdfcairo
set output "Punto1A.pdf"
set xrange [0.0:80.0]
set yrange [0.0:1.5]
set xlabel "Time"
set ylabel "SIR"
set title "SIR Model with beta = 0.35, gamma = 0.08"
plot "Punto1A.txt" u 2:3 w l lw 2.5 t "Susceptibles", "Punto1A.txt" u 2:4 w l lw 2.5 t "Infected", "Punto1A.txt" u 2:5 w l lw 2.5 t "Recovered"