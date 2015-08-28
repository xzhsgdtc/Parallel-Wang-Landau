#! /usr/bin/env python

import copy
import MyMath
import numpy as np
from sys import stdout
from sys import stdin
from sys import stderr
from sys import exit


fmt_1d_histogram = "%8d   %8.5f   %8.5f   %10d    %10.5f"
fmt_2d_histogram = "%8d   %8.5f   %8.5f   %8.5f   %8.5f   %10d    %10.5f"

class HistogramND:
    def __init__(self, dim):
        self.dim                = dim
        self.num_bins           = 0
        self.bins               = None
        self.flag               = np.zeros(dim, dtype=np.int)
        self.num_bins_each_dim  = np.zeros(dim, dtype=np.int)
        self.weight_each_dim    = np.zeros(dim)
        self.step_each_dim      = np.zeros(dim)
        self.max_each_dim       = np.zeros(dim)
        self.min_each_dim       = np.zeros(dim)

    def copy(self, his):
        self.dim                = copy.deepcopy(his.dim)
        self.num_bins           = copy.deepcopy(his.num_bins)
        self.bins               = copy.deepcopy(his.bins)
        self.flag               = copy.deepcopy(his.flag)
        self.num_bins_each_dim  = copy.deepcopy(his.num_bins_each_dim)
        self.weight_each_dim    = copy.deepcopy(his.weight_each_dim)
        self.step_each_dim      = copy.deepcopy(his.step_each_dim)
        self.max_each_dim       = copy.deepcopy(his.max_each_dim)
        self.min_each_dim       = copy.deepcopy(his.min_each_dim)
        return self

    def dos(self):
        if self.dim == 1:
            one = np.append(np.reshape(self.bins[:,1],(self.num_bins,1)),np.reshape(self.bins[:,2*self.dim+2],(self.num_bins,1)),axis=1)
        elif self.dim == 2:
            one = np.append(np.reshape(self.bins[:,1],(self.num_bins,1)),np.reshape(self.bins[:,3],(self.num_bins,1)),axis=1)
            one = np.append(one, np.reshape(self.bins[:,2*self.dim+2],(self.num_bins,1)),axis=1)
        return one

    def norm_dos(self):
        if self.dim == 1:
            one = np.append(np.reshape(self.bins[:,1],(self.num_bins,1)),np.reshape(self.bins[:,2*self.dim+2],(self.num_bins,1)),axis=1)
            one -= (0, one[0,1])
        elif self.dim == 2:
            one = np.append(np.reshape(self.bins[:,1],(self.num_bins,1)),np.reshape(self.bins[:,3],(self.num_bins,1)),axis=1)
            one = np.append(one, np.reshape(self.bins[:,2*self.dim+2],(self.num_bins,1)),axis=1)
            one -= (0, 0,one[0,2])
        return one

    def read(self, filename=None):
        #read in parameters
        with open(filename, 'r') as fin:
            for line in fin:
                line = line.strip()
                words = line.split(':')
                if len(words) == 0:
                    continue
                if words[0] == '# Dimension':
                    dim = int(words[1])
                    if self.dim != dim:
                        stderr.write("\n*** Error: dimension %d != %d !\n\n" % (dim,self.dim))
                        exit(1)
                elif words[0] == '# NumberofBins':
                    self.num_bins = int(words[1])
                elif words[0] == '# NumberOfBinsEachDimension':
                    w = words[1].split()
                    for i in xrange(0,self.dim):
                        self.num_bins_each_dim[i] = w[i]
                elif words[0] == '# WeightOfEachDimension':
                    w = words[1].split()
                    for i in xrange(0,self.dim):
                        self.weight_each_dim[i] = w[i]
                elif words[0] == '# StepOfEachDimension':
                    w = words[1].split()
                    for i in xrange(0,self.dim):
                        self.step_each_dim[i] = w[i]
                elif words[0] == '# MaximumOfEachDimension':
                    w = words[1].split()
                    for i in xrange(0,self.dim):
                        self.max_each_dim[i] = w[i]
                elif words[0] == '# MinimumOfEachDimension':
                    w = words[1].split()
                    for i in xrange(0,self.dim):
                        self.min_each_dim[i] = w[i]

        self.bins = np.loadtxt(filename)
        self.num_bins = self.bins.shape[0]

        for i in xrange(0,self.dim):
            if np.allclose( self.bins[:,i*2+1], self.bins[:,i*2+2] ):
               self.flag[i] = 1
            else:
               self.flag[i] = 0

    def write(self, filename=None):
        if filename == None:
            stderr.write("\n***Please specify the filename to write!\n\n")
            exit(1)
        if self.dim == 1:
            np.savetxt(filename, self.bins, fmt=fmt_1d_histogram)
        elif self.dim == 2:
            np.savetxt(filename, self.bins, fmt=fmt_2d_histogram)

        with open(filename,'a') as fout:
            fout.write("\n# Dimension: %d \n" % self.dim)
            fout.write("\n# NumberofBins: %d \n" % self.num_bins)
            temp = ""
            for i in xrange(0, self.dim-1):
                temp += ' %d ' % self.num_bins_each_dim[i]
            temp += ' %d \n' % self.num_bins_each_dim[self.dim-1]
            fout.write("\n# NumberOfBinsEachDimension: " + temp)

            temp = ""
            for i in xrange(0, self.dim-1):
                temp += ' %d ' % self.weight_each_dim[i]
            temp += ' %d \n' % self.weight_each_dim[self.dim-1]
            fout.write("\n# WeightOfEachDimension:" + temp)

            temp = ""
            for i in xrange(0, self.dim-1):
                temp += ' %d ' % self.step_each_dim[i]
            temp += ' %d \n' % self.step_each_dim[self.dim-1]
            fout.write("\n# StepOfEachDimension: " + temp)

            temp = ""
            for i in xrange(0, self.dim-1):
                temp += ' %d ' % self.max_each_dim[i]
            temp += ' %d \n' % self.max_each_dim[self.dim-1]
            fout.write("\n# MaximumOfEachDimension: " + temp)

            temp = ""
            for i in xrange(0, self.dim-1):
                temp += ' %d ' % self.min_each_dim[i]
            temp += ' %d \n' % self.min_each_dim[self.dim-1]
            fout.write("\n# MinimumOfEachDimension: " + temp)

    # split 2D histogram into pieces of 1D histogram
    # @param split_dim  split according to which dimension
    def split(self):
        split_dim = 1
        if self.dim == 1:
            return [copy.deepcopy(self)]
        num = self.num_bins_each_dim[split_dim]
        his_1d = []
        for i in xrange(int(num)):
            his = HistogramND(2)
            for d in xrange(self.dim):
                his.flag[d]                =    self.flag[d]
                his.num_bins_each_dim[d]   =    self.num_bins_each_dim[d]
                his.weight_each_dim[d]     =    self.weight_each_dim[d]
                his.step_each_dim[d]       =    self.step_each_dim[d]
                his.max_each_dim[d]        =    self.max_each_dim[d]
                his.min_each_dim[d]        =    self.min_each_dim[d]
            his.num_bins =    self.num_bins_each_dim[(split_dim+1)%2]
            his.num_bins_each_dim[split_dim] = 1
            his.bins = copy.deepcopy(self.bins[i*his.num_bins : (i+1)*his.num_bins, : ])
            his.max_each_dim[split_dim] = his.bins[0,split_dim*2+2]
            his.min_each_dim[split_dim] = his.bins[his.num_bins-1,split_dim*2+1]
            his_1d.append(copy.deepcopy(his))
        return his_1d

    # merge two histogram together according to their first derivatives
    # @param h histogram instance
    def merge(self, h):
        if self.min_each_dim[0] < h.min_each_dim[0] and self.max_each_dim[0] < h.max_each_dim[0]    \
                and ( self.dim ==1 or self.dim == 2 ):
            one = np.append(np.reshape(self.bins[:,1],(self.num_bins,1)),np.reshape(self.bins[:,2*self.dim+2],(self.num_bins,1)),axis=1)
            two = np.append(np.reshape(h.bins[:,1],(h.num_bins,1)),np.reshape(h.bins[:,2*h.dim+2],(h.num_bins,1)),axis=1)
        else:
            stderr.write("\n*** Error: range of histogram is invalid! Merge operation is illegal!\n\n")
            exit(1)

        deriv_one = MyMath.derivatives(one)
        deriv_two = MyMath.derivatives(two)

        start_index = -1
        for i in xrange(0,self.num_bins-4):
            if np.allclose(deriv_one[i,0],deriv_two[0,0]):
                start_index = i
                break
        end_index = -1
        for i in xrange(0,h.num_bins-4):
            if np.allclose(deriv_one[self.num_bins-5, 0],deriv_two[i, 0]):
                end_index = i
                break

        if start_index == -1 or end_index == -1:
            stderr.write("\n*** Error: histograms are not compatible! Merge Operation is illegal!\n\n")
            exit(1)

        deriv_one = deriv_one[start_index :           , :]
        deriv_two = deriv_two[            :end_index+1 , :]

        diff = np.abs(deriv_one[:,2]-deriv_two[:,2])
        min_pos = np.argmin(diff)
        min_pos_in_one = min_pos + 2 + start_index
        min_pos_in_two = min_pos + 2


        temp = np.copy(h.bins)
        if self.dim == 1:
            offset = self.bins[min_pos_in_one, 4]
            temp += (0,0,0,0, offset - temp[min_pos_in_two,4])
        elif self.dim == 2:
            offset = self.bins[min_pos_in_one, 6]
            temp += (0,0,0,0,0,0, offset - temp[min_pos_in_two,6])

        self.bins = np.append(self.bins[:min_pos_in_one,:], temp[min_pos_in_two:,:], axis=0)
        self.num_bins = self.bins.shape[0]
        self.num_bins_each_dim[0] = self.bins.shape[0]
        self.max_each_dim = np.copy(h.max_each_dim)

        for i in range(0,self.bins.shape[0]):
            self.bins[i,0] = i

    def __str__(self):
        print ""
        print "# Dimension   %d " % self.dim
        temp = " Each Dim  "
        for i in xrange(0, self.dim-1):
            temp += ' %d, ' % self.num_bins_each_dim[i]
        temp += ' %d ' % self.num_bins_each_dim[self.dim-1]
        print "# Number of Bins  %d " % self.num_bins + temp
        temp = ""
        for i in xrange(0, self.dim-1):
            temp += ' %d, ' % self.weight_each_dim[i]
        temp += ' %d ' % self.weight_each_dim[self.dim-1]
        print "# Weight Each Dim " + temp
        temp = ""
        for i in xrange(0, self.dim-1):
            temp += ' %d, ' % self.step_each_dim[i]
        temp += ' %d ' % self.step_each_dim[self.dim-1]
        print "# Step Each Dim " + temp
        temp = ""
        for i in xrange(0, self.dim-1):
            temp += ' %d, ' % self.max_each_dim[i]
        temp += ' %d ' % self.max_each_dim[self.dim-1]
        print "# Maximum Each Dim " + temp
        temp = ""
        for i in xrange(0, self.dim-1):
            temp += ' %d, ' % self.min_each_dim[i]
        temp += ' %d ' % self.min_each_dim[self.dim-1]
        print "# Minimum Each Dim " + temp
        print ""
        print self.bins
        return ''

if __name__=='__main__':

    from sys import argv
    h = HistogramND(2)
    filename = argv[1]
    h.read(filename)
    print h.dos()
    #from sys import argv
    #import os
    #import DOS

    #if len(argv) != 5:
    #    print """
    #    %s [folder path] [file format] [start index] [end index]
    #    """ % argv[0]
    #    exit(1)

    #path = argv[1]
    #fmt = argv[2]
    #index = [int(argv[3]), int(argv[4])]

    #num = index[1] - index[0] + 1
    #print "# Number of file    %d" % num

    #hs = []
    #for i in range(index[0],index[1]+1):
    #    h = HistogramND(2)
    #    filename = os.path.join(path,fmt%i)
    #    h.read(filename)
    #    h = h.split()
    #    hs.append(h)

    #h = hs[0][0]
    #for i in range(1,num):
    #    h.merge(hs[i][0])
    #h.write("temp.test")
    #a = h.dos()
    #a = a[35:,:]

    #dT = 0.05
    #num_points=100;
    #thermo_result = DOS.thermo(a, dT, num_points, 1000)

    #np.savetxt("thermo.test.dat",thermo_result)
