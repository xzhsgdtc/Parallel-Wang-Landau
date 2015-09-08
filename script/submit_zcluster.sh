#!/bin/bash

# submit script for submitting parallel jobs to zcluster
# this script should be submitted to zcluster through following command
# qsub -q rcc-30d -pe mpi [number of processes] [script filename]
# e.g. qsub -q rcc-30d -pe mpi 5 submit_zcluster.sh

export LD_LIBRARY_PATH=/usr/local/mpich2/1.4.1p1/pgi123/lib:${LD_LIBRARY_PATH}

# please customize below submit command if needed
# time mpirun -np $NSLOTS ./PWLSIM  [input file]  [dimension of DOS]  [ ... wins for each dim ... ] [ ... overlap for dim ... ]
time mpirun -np $NSLOTS ./PWLSIM wang_landau.input 1 5 0.75 

