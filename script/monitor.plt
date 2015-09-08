#! /usr/bin/gnuplot


ifname="../record/P0000000_Track_0001.dat"


set style line 1 lc rgb "red" pt 6
set style line 2 lc rgb "black" pt 7


set multiplot layout 2,1

#set yrange[-5000:-3000]
set xtics nomirror
set ytics nomirror
plot ifname u 1:6 w l notitle  ls 1
#set yrange[94:97]
plot ifname u 1:7 w lp notitle  ls 2

unset multiplot

pause 2
reread
