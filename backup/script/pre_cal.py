#! /usr/bin/env python


from sys import argv, exit

if len(argv) != 4:
    print """
        Usage: %s [total bins] [No. of bins each process] [overlap]
    """ % argv[0]
    exit()


num_bins = int(argv[1])
overlap = float(argv[3])
bins_each = int(argv[2])

num_procs = float(1 + (num_bins - bins_each) / ((1-overlap) * bins_each))
overlap_bins = float(bins_each*overlap)

print """
    Input:

        *Number of Total Bins = %d
        *Number of Bins Each Proc. = %d
        *Overlap Fraction = %f

    Output:

        *Number of process = %f
        *Overlap Bins = %f
""" % (num_bins, bins_each, overlap, num_procs, overlap_bins)
