set terminal qt size 800,800

set title "Pseudo Cell Movement on Curved Surface"
set xlabel "x"
set ylabel "y"
set zlabel "z"

set xrange [-4:4]
set yrange [-4:4]
set zrange [-0.4:0.4]

set grid
set border 4095
unset key

set xyplane at 0
set view 70, 225

set samples 80
set isosamples 80

surface(x,y) = 0.2*cos(x)

# Gradient surface
set pm3d at s
set palette rgbformulae 22,13,-31
set style fill transparent solid 0.35

Nsteps = 10

do for [i=0:Nsteps-1] {
    set title sprintf("Pseudo Cell Movement on Curved Surface - Frame %d", i)

    splot \
        surface(x,y) with pm3d, \
        'membrane.dat' index i using 1:2:($3 + 0.02) with points pt 7 ps 1.4 lc rgb "black", \
        'interior.dat' index i using 1:2:($3 + 0.02) with points pt 7 ps 1.0 lc rgb "red"
        
    pause .5
}
