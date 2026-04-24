set terminal wxt size 800,800
set title "Pseudo Cell Movement in 3D"
set xlabel "x"
set ylabel "y"
set zlabel "z"

set xrange [-4:4]
set yrange [-4:4]
set zrange [-1:1]

set xyplane at 0
set view 60,30
set view equal xy

set size 1,1
set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0

set grid
set border 4095
unset key

Nsteps = 10

do for [i=0:Nsteps-1] {
    splot \
        'membrane.dat' index i using 1:2:3 with points pt 7 ps 1.0 lc rgb "black", \
        'interior.dat' index i using 1:2:3 with points pt 7 ps 0.6 lc rgb "red"

    pause .5
}