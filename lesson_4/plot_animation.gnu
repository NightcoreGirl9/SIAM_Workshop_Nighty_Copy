set terminal gif animate delay 50 size 800,800
set output "pseudo_movement.gif"

set title "Pseudo Cell Movement"
set xlabel "x"
set ylabel "y"

set xrange [-4:4]
set yrange [-4:4]

set size square
set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0

set grid
set border 4095
unset key

Nsteps = 10

do for [i=0:Nsteps-1] {
    set title sprintf("Pseudo Cell Movement - Frame %d", i)

    plot \
        'membrane.dat' index i using 1:2 with points pt 7 ps 1.0 lc rgb "black", \
        'interior.dat' index i using 1:2 with points pt 7 ps 0.6 lc rgb "red"
}

set output