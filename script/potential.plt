#! /usr/bin/gnuplot -persist

set term x11

t = 1**(1/6)

step(x)= x > t ? 0 : 1 

set xrange[0.7:1.5]
plot -0.5*30*1.3*1.3*log(1-(x/1.3)**2) + 4*1*((1/x)**12 - (1/x)**6) + step(x) ,\
     -0.5*40*0.09*log(1-((x-1)/0.3)**2)
#plot 4*1*((1/x)**12 - (1/x)**6)
