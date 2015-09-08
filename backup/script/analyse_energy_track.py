#! /usr/bin/env python

import numpy as np
from sys import argv, exit
import matplotlib.pyplot as plt




if __name__ == "__main__":
    if len(argv) != 8:
        print """
        usage: %s [energy tracking file] [energy min] [energy high] [bin width] [lipid min] [lipid max] [lipids]
        """ % argv[0]
        exit(1)

    ifile = argv[1]
    energy_range = [float(argv[2]), float(argv[3])]
    lipids_range = [int(argv[5]), int(argv[6])]
    bin_width = float(argv[4])
    num_lipid = int(argv[7])

    print """
    ######################################

    input file:  %s
    energy range: %f,   %f
    lipids range: %d,   %d
    bin width: %f

    """ % (ifile, energy_range[0], energy_range[1], lipids_range[0], lipids_range[1], bin_width)

    content = np.loadtxt(ifile)
    num_lines = content.shape[0]
    content = np.reshape(content[:,5:7],(num_lines,2))

    data = np.zeros((1,2))
    for i in range(num_lines):
        if content[i,1] == num_lipid:
            data = np.append(data,np.reshape(content[i,:], (1,2)),axis=0)

    his = np.histogram(data, bins=100, range=energy_range)
    np.savetxt("test.dat",his[0], "%d")
    def show_his(data):
        plt.hist(data[:,0],100,range=energy_range)
        plt.show()

    show_his(data)


