#! /usr/bin/env python

import math
from sys import argv, exit

if len(argv) != 4 :
    print """
    usage: %s [alpha] [start T] [threshold]
    """
    exit(1)

alpha = float(argv[1])
start_t = float(argv[2])
stop_t = float(argv[3])

T = start_t
count = 0
while T > stop_t:
    T *= alpha
    count +=1

print "Total Count = %d" %count
