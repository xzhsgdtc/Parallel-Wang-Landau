#!/bin/sh

### username
USER=sgjerry

### Tell queue system to use the directory of the submit script ###
### as the current working directory ###
#$ -cwd

### Tell queue system to submit this job to a specific queue ###
###(currently imac, opteron, or xeon ) ###
#$ -q xeon
#$ -pe mpi 12

### make results directory in your home on hal ###
### mkdir $JOB_ID
resultdir=`pwd`

### make directory on local machine in scratch ###
#mkdir -p /scratch/$USER/$JOB_ID

### Tell queue system to write it's output and error files to the scratch directory ###
### PLEASE NOTE since these output files are created before this script is run, ###
### they must be in /scratch/USERname since /scratch/username/$JOB_ID doesn't exist yet ###
#$ -o ./$JOB_NAME.o$JOB_ID

#$ -e ./$JOB_NAME.e$JOB_ID

### execute code ###
time mpirun -np 12 ./PWLSIM wang_landau.input 1 12 0.8 0 928749
 

###copy files back to home directory located on hal ###
### cp -r * $resultdir

###clean up scratch directory since space is limited on local machines ###
#rm -rf /scratch/sgjerry/$JOB_ID

