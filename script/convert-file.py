
####################################################
#                                                  #
#                                                  #
#           This script is designed to convert     #
#           between different file formats         #
#           @Author Jerry Shi                      #
#           @Version 1.0 June 15, 2014             #
#               mol2 -> pdb                        #
#                                                  #
#                                                  #
####################################################

#! /usr/bin/env python

import math
from mylib.model import GeneralModel
from mylib.filehelper import reader


def fun_mol2_to_pdb(ifile, ofile, shift=None):

    models = reader.read_mol2(ifile)
    num = len(models)

    print "    -> number of models:    %d"  % num
    for i in xrange(0, num):
        m = models[i]

        if shift !=None :
            m.shift_and_apply_periodic(0, shift[1], shift[0])
            m.shift_and_apply_periodic(1, shift[2], shift[0])
            m.shift_and_apply_periodic(2, shift[3], shift[0])
        m.write(ofname, 'pdb', i)


if __name__ == '__main__':
    import argparse
    from sys import argv, exit
    from mylib.filehelper import pdb

    parser = argparse.ArgumentParser(description="This script is designed to convert file fromat.")
    parser.add_argument('ifile',    nargs=1, type=str,   help="input filename")
    parser.add_argument('ofile',    nargs=1, type=str,   help="output filename")
    parser.add_argument('format',   nargs=2, type=str,   help = "[from format] [to format]")
    parser.add_argument('-shift',   nargs=4, type=float, help ="box-length, dim-0-shift, dim-1-shift, dim-2-shift")

    args = parser.parse_args()

    ifname = args.ifile[0]
    ofname = args.ofile[0]
    from_format = args.format[0]
    to_format = args.format[1]
    shift = args.shift


    print """
    ########################## parameters ##########################
    * input-file: %s
    * output-file: %s
    * file-format: from %s  to %s
    """ % (ifname, ofname, from_format, to_format)


    if from_format == 'mol2' and to_format == 'pdb':
        fun_mol2_to_pdb(ifname, ofname, shift);
    else:
        print "Error! file formats are not supported!"
        exit(1);


