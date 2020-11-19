#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1 
#PBS -N install
#PBS -j oe  
#PBS -m ae   
#PBS -A PAS1627

#
# The following lines set up the R environment
#
#module load intel
module load openmpi
#module load mkl
module load R/4.0.2-gnu9.1
#
# Move to the directory where the job was submitted
#
cd $PBS_O_WORKDIR
# parallel R: submit job with one MPI master
mpirun -np 1 R --slave < 00_preparation.R
