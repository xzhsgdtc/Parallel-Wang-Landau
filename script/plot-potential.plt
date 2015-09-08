#! /usr/bin/gnuplot -persist

set term x11


lj_sigma = 0.890899
lj_epsilon = 1.000000
lj_cutoff = 2.227247
lj_shift = 4*lj_epsilon*(lj_sigma**12/(lj_cutoff**12) - lj_sigma**6/(lj_cutoff**6))
#
#set xrange[0.5:2.227247]
#plot 4*lj_epsilon*(lj_sigma**12/(x**12) - lj_sigma**6/(x**6))- lj_shift

fene_k  = 40.000000
fene_r  = 0.300000
fene_ro = 1.000000

set xrange [0.6:1.4]
wca (r) = r > 1 ? 0 :  (4*lj_epsilon*(lj_sigma**12/(r**12) - lj_sigma**6/(r**6)) + lj_epsilon)
fene(x) =  -fene_k * (fene_r**2)*0.5*log((1-((x-fene_ro)/fene_r)**2))  
print fene(10 )
