#! /usr/bin/env python

import numpy as np

fene_k  = 40.000000
fene_r  = 0.300000
fene_ro = 1.000000

def fene(x):
    return -fene_k * (fene_r**2)*0.5*np.log((1-((x-fene_ro)/fene_r)**2))


for i in xrange(1,100):
    x = 0.6 + i*0.01
    print x, fene(x)
