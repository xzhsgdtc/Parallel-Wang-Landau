#! /usr/bin/gnuplot -persist


aa=0.2
bb=0.5

Popt=0.36


set xrange[0:1]

plot log(aa*Popt+bb)/log(aa*x+bb)
