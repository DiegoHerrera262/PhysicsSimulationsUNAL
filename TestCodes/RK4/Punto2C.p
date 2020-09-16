set terminal pdfcairo
set output "SomeDrumEigenfunc.pdf"
set ytics 1
set xtics 0.2
set grid
set xrange [0.0:1.0]
set yrange [-0.5:1.3]
set xlabel "r"
set ylabel "R(r)"
set title "Drum Eigenfunctions - Count Nodes to Identify"
plot "DrumEigen1.txt" u 2:3 w l lw 2.5 t "1st", "DrumEigen2.txt" u 2:3 w l lw 2.5 t "2nd","DrumEigen3.txt" u 2:3 w l lw 2.5 t "3rd" ,"DrumEigen4.txt" u 2:3 w l lw 2.5 t "4th" ,"DrumEigen5.txt" u 2:3 w l lw 2.5 t "5th"
