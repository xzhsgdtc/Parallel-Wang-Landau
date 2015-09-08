#! /usr/bin/env python

from ..model import GeneralModel

def read_mol2(ifile):
    num_atoms = -1
    num_bonds = -1
    read_status = "";
    models = []
    with open(ifile, 'r') as fin:
        content = fin.readlines()
        size = len(content)
        i = 0;
        while i < size:
            words = content[i].split()
            if len(words) == 0:
                i +=1;
                continue

            if words[0] == '@<TRIPOS>MOLECULE':
                line = content[i+2]
                words = line.split();
                num_atoms = int(words[0])
                num_bonds = int(words[1])
                m = GeneralModel(num_atoms, num_bonds)
                i += 5
            elif words[0] == '@<TRIPOS>ATOM':
                for j in xrange(0, num_atoms):
                    line = content[i+j+1]
                    words = line.split()
                    if j+1 != int(words[0]):
                        print "Error in reading atom!"
                        exit(1)
                    m.add_atom(j+1, words[1], words[2:5])

                i += num_atoms
            elif words[0] == '@<TRIPOS>BOND':
                for j in xrange(0, num_bonds):
                    line = content[i+j+1]
                    words = line.split()
                    if j+1 != int(words[0]):
                        print "Error in reading bond!"
                        exit(1)
                    m.add_bond(words[0], words[1], words[2])
                i += num_bonds
                models.append(m)
            else:
                i += 1
    return models
