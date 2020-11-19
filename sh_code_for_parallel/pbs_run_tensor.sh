#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1
#PBS -l mem=40gb	
#PBS -N spibayes_tensor
#PBS -j oe  
#PBS -m ae   
#PBS -A PAS1627

#
# The following lines set up the R environment
#
#module load intel
module load openmpi
#module load mkl
module load netcdf
module load R/4.0.2-gnu9.1
#
# Move to the directory where the job was submitted
#
cd $PBS_O_WORKDIR
# parallel R: submit job with one MPI master
mpirun -np 1 R --slave < 10-gauge_tensor_output.R
