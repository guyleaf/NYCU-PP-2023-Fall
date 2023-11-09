#!/usr/bin/gnuplot

view_type = "view2"
data = "q1_speed_up_".view_type.".dat"
image = "q1_speed_up_".view_type.".png"

set title view_type." - Speed Up, multi-threaded vs single-threaded"
set xlabel "thread"
set ylabel "ratio"
set xtics 1

set grid

# visualize on vscode-previewer
plot data using 1:2 with lines title ""

# export image
set term png enhanced font "Verdana,10"
set output image
replot