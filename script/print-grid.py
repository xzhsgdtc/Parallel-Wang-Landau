#! /usr/bin/env python

from sys import stdout

num_energy_wins = 21
num_lipids_wins = 36


for i in range(num_lipids_wins):
    for j in range(num_energy_wins):
        stdout.write( "%4d " % ( i*num_energy_wins + j))
    print
