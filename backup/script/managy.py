#! /usr/bin/env python

from sys import argv, exit
import random
import numpy as np
import os, shutil


dbpath ="../../../ResearchCSP/OLModel/Data_N_Analysis/RawData/Stampede/ConfsDatabase"
opath = ""
nfformat = "P%07d_InitConf.mol2"


if __name__=="__main__":
    if len(argv) != 9:
        print """
        usage: %s [EMin.] [EMax.] [NMin.] [NMax] [No. EWins.] [No. NWins.] [Eoverlap] [Noverlap]
        Note, all are inclusive
        """ % argv[0]
        exit(1)

    Emin = float(argv[1])
    Emax = float(argv[2])
    Nmin = int(argv[3])
    Nmax = int(argv[4])
    NoEWins = int(argv[5])
    NoNWins = int(argv[6])
    Eoverlap = float(argv[7])
    Noverlap = float(argv[8])
    total_processes = NoEWins*NoNWins
    Erange = Emax - Emin
    Nrange = Nmax - Nmin + 1
    Ewidth = Erange / (1 + (1-Eoverlap)*(NoEWins-1))
    Nwidth = Nrange / (1 + (1-Noverlap)*(NoNWins-1))

    print"""
      #Energy        %6.3f  to  %6.3f
      #Lipid No.       %5d    to  %5d
      #Win No.       E %5d    , N %5d

      #Total processes   %d
      #Erange -- Ewidth    %6.3f  --  %6.3f
      #Nrange -- Nwidth        %5d    --  %5d
    """ % (Emin, Emax, Nmin, Nmax, NoEWins, NoNWins, total_processes, Erange, Ewidth, Nrange, Nwidth)


    for ni in range(0, NoNWins):
        for ei in range(0, NoEWins):
            myEmin = Emin + (1-Eoverlap)*ei*Ewidth
            myEmax = myEmin + Ewidth
            myNmin = Nmin + (1-Noverlap)*ni*Nwidth
            myNmax = myNmin + Nwidth
            pid = ni*NoEWins+ei
            print "process %4d  Min %6.3f %5d  Max %6.3f %5d" % (pid, myEmin, myNmin, myEmax, myNmax)
            E = int(np.abs(random.random()*(myEmax-myEmin) + myEmin))
            N = int(random.random()*(myNmax-myNmin) + myNmin)
            folderpath = os.path.join(dbpath, "%d_%d" % (E,N))
            while not os.path.exists(folderpath):
                E = int(np.abs(random.random()*(myEmax-myEmin) + myEmin))
                N = int(random.random()*(myNmax-myNmin) + myNmin)
                folderpath = os.path.join(dbpath, "%d_%d" % (E,N))

            print "Randomly pick Energy = %6d, Lipids = %5d" % (E,N)
            nfname = nfformat % pid
            fs = os.listdir(folderpath)
            shutil.copyfile(fs[random.randint(0, len(fs)-1)], os.path.join(opath, nfname))


