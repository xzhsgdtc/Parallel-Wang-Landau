#! /usr/bin/python2.7


##################################################
###
#
#   This script is desinged to track the position
#   of each package from data files of Parallel WL
#   simulation
#
#   @Author Jerry Shi
#   @Data   May 29th
#
###
##################################################


import os, sys
import numpy as np


if len(sys.argv) != 6:
    print """
    Usage: ./script     [Input file format: e.g. P0000%02d.energy_track.dat]    [No. of Packages (i.e. No. of Procs)]    [Start Index]  [No. of Wins in X] [No. of Wins in Y]
    """
    sys.exit(1)


fileformat = sys.argv[1]
num_packages = int(sys.argv[2])
start_index = int(sys.argv[3])
num_wins = (int(sys.argv[4]), int(sys.argv[5]))
num_procs_win = num_packages/(num_wins[0]*num_wins[1])

print """
######################################
### Processing the track of package

    Input file format: %s
    No. of Input file: %d (from %d)
    No. of Procs Per Window: %d
    No. of windows (X,Y): (%d, %d)
    _________________________________________
""" % (fileformat, num_packages, start_index, num_procs_win, num_wins[0], num_wins[1])


files=[]
tracks=[]
#counter= np.zeros((num_packages,1))
print """    readinng data ..."""


max_lines=-1
min_lines=-1
for i in range(start_index, start_index+num_packages):
    filename = fileformat % i
    content = np.loadtxt(filename)
    lines = content.shape[0]

    if max_lines == -1 or max_lines < lines:
        max_lines = lines
    if min_lines == -1 or min_lines > lines:
        min_lines = lines

    content = np.reshape(content[:,2:4],(lines,2))
    content = content.astype(np.int64)
    files.append(np.copy(content))

for i in range(start_index, start_index+num_packages):
    tracks.append(-np.ones((min_lines,content.shape[1]+3)))      # line number, package id, package position, pos-x, pos-y
print "     Maximum lines: %d" % max_lines
print "     Minimum lines: %d" % min_lines


j=0
while j < min_lines:
    for i in range(start_index, start_index+num_packages):
        if j < files[i].shape[0]:
            pid = int(files[i][j,0])
            pos = int(files[i][j,1])
            pline = pos / (num_wins[0]*num_procs_win)
            prow = (pos/num_procs_win) % num_wins[0]
            tracks[pid][j] = [j, pid, pos, pline, prow]
     #       counter[pid] = j
    j += 1

print """    Output file ..."""

for i in range(start_index, start_index+num_packages):
    filename = "P%07d.package.track" % i
    #tracks[i] = tracks[i][0:int(counter[i]+1),:]
    np.savetxt(filename,tracks[i],fmt="%d    %d    %d     %d      %d")

print """
DONE !
#######################################
"""
