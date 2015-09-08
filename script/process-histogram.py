#! /usr/bin/env python

import  os
from sys import path, argv, exit
path.append("/csphome/sgjerry/git/parallel-lipid-model/script/Library")
import DOS
from HistogramND import HistogramND
import argparse, random
import numpy as np
import MyMath

def split_files():

    path = argv[2]
    fmt = argv[3]
    index = [int(argv[4]), int(argv[5])]

    num = index[1] - index[0] + 1
    print "# Number of file    %d" % num


    for i in range(index[0],index[1]+1):
        h = HistogramND(2)
        filename = os.path.join(path,fmt%i)
        print "    -> %s" % filename
        h.read(filename)
        h = h.split()
        for one in h:
            fn = "P%07d_%d.dat" % (i, one.min_each_dim[1])
            fn = os.path.join(path,fn)
            one.write(fn)

def resample(hs_files):
    f = []
    for one in hs_files:
        i = (int(random.random()*1000000)) % len(one)
        f.append(one[i])
    return f

def merge_files():

    path = argv[2]
    fmt = argv[3]
    wins = [int(argv[4]), int(argv[5])]
    lipids = int(argv[6])
    num_resamples = int(argv[7])

    hs_files = []
    for i in range(0, wins[0]):
        hs_files.append([])
        for j in range(0, wins[1]):
            filename = fmt % ((i+j*wins[0]), lipids)
            filename = os.path.join(path, filename)
            if os.path.exists(filename):
                hs_files[i].append(filename)

    result = None
    for r in range(0,num_resamples):
        h = HistogramND(2)
        fs = resample(hs_files)
        h.read(fs[0])
        for i in range(1,len(fs)):
            th = HistogramND(2)
            th.read(fs[i])
            h.merge(th)
        h = h.norm_dos()
        if result == None:
            result = np.copy(  np.append( np.reshape(h[:,0], (h.shape[0],1)), np.reshape(h[:,2], (h.shape[0],1)) , axis=1))
        else:
            result = np.append(result, np.reshape(h[:,2],(result.shape[0],1)), axis=1)

    opath = "."
    final = np.reshape(result[:,0],(result.shape[0],1))
    final = np.append(final,np.reshape(np.mean(result[:,1:], axis=1),(result.shape[0],1)), axis=1)
    final = np.append(final,np.reshape(np.std(result[:,1:], axis=1), (result.shape[0],1)),axis=1)
    np.savetxt(os.path.join(opath,"Lipid%d.dos_deriv.dat" % lipids), MyMath.derivatives(final[:,0:2]) )

    array_lipid = np.ones((final.shape[0],1)) * lipids
    final = np.append(final,np.reshape(array_lipid,((final.shape[0],1))),axis=1)
    np.savetxt(os.path.join(opath,"Lipid%d.ave_dos.dat" % lipids), final, fmt="%.2f   %.2f   %.2f  %d")

    dT = 0.01
    kb=1
    num_files = 1;
    num_points=1000;
    num_atoms=1000

    result = result[67:, ]

    thermo_result = DOS.thermo(result, dT, num_points, 1000)

    np.savetxt(os.path.join(opath,"Lipid%d.thermo_dos.dat" % lipids),thermo_result)


if __name__ == "__main__":
    from sys import argv, exit, stderr
    if len(argv) < 2:
        print """
        %s [function: split or merge] [ ... ]
        """ % argv[0]
        exit(1)

    func = argv[1]
    if func == 'split':
        if len(argv) != 6:
            print """
        split function: %s split [folder path] [file format] [start index] [end index]
            """ % argv[0]
            exit(1)
        split_files()
    elif func == 'merge':
        if len(argv) != 8:
            print """
        merge function: %s merge [folder path] [file format] [Energy Dim. Wins] [Lipids Dim. Wins] [lipid number] [number of resamples]
            """ % argv[0]
            exit(1)
        merge_files()
