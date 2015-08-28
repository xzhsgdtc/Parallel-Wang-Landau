#! /usr/bin/env python


import numpy as np
import math
from filehelper import pdb

class GeneralModel:

    def __init__(self, num_atoms, num_bonds):
        self.n_atoms = num_atoms
        self.n_bonds = num_bonds
        self.atoms_coord = np.zeros((self.n_atoms, 3))
        self.atoms_type  = np.chararray(self.n_atoms)
        self.bonds       = np.zeros( (self.n_bonds,2), dtype=np.int)

    def add_atom(self, i, t, coord):
        i = int(i)
        self.atoms_coord[i-1] = np.copy(np.asarray(coord))
        self.atoms_type[i-1] = t

    def add_bond(self, i, a1, a2):
        i = int(i)
        self.bonds[i-1] = np.asarray([a1, a2])

    def shift_and_apply_periodic(self, d, shift, box_length):
        for i in xrange(0,self.n_atoms):
            self.atoms_coord[i][d] -= shift
            self.atoms_coord[i][d] = self.atoms_coord[i][d] - box_length * math.floor(self.atoms_coord[i][d]/box_length)

    def write(self, ofname, oftype, mi=1):
        if oftype == 'pdb':
            with open(ofname, 'aw') as fout:
                fout.write(pdb.beg_model(mi) + "\n")
                for i in xrange(0, self.n_atoms):
                    wline = pdb.atom (      i+1,  self.atoms_type[i], '', self.atoms_type[i], 'A', mi, '', \
                                            self.atoms_coord[i][0], self.atoms_coord[i][1], self.atoms_coord[i][2],\
                                            1.00, 0.00, self.atoms_type[i], '' )
                    fout.write(wline+"\n")
                fout.write(pdb.ter(self.n_atoms +1, 'E', 'A',1, '') + '\n')
                for i in xrange(0, self.n_bonds):
                    wline = pdb.connect(self.bonds[i][0], self.bonds[i][1])
                    fout.write(wline+"\n")
                fout.write(pdb.end_model() + "\n")

    def __distance(self, i, j):
        d = self.atoms_coord[i] - self.atoms_coord[j]
        return np.sqrt(np.sum(np.square(d)))

    def bond_distance(self):
        print "#Number of bonds:   %d" % self.n_bonds
        for i in xrange(0, self.n_bonds):
            print "Bond  %4d     Monomer %4d -- Monomer %4d   ->  %f" % (i+1, self.bonds[i][0], self.bonds[i][1], self.__distance(self.bonds[i][0]-1,self. bonds[i][1]-1))


