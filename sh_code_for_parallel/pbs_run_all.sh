#!/bin/bash
#PBS -l walltime=15:00:00
#PBS -l nodes=1:ppn=1 
#PBS -N spibayes_paper
#PBS -j oe  
#PBS -m ae   
#PBS -A PAS1627

#
# The following lines set up the R environment
#
module load intel
module load openmpi
module load mkl
module load R
#
# Move to the directory where the job was submitted
#
cd $PBS_O_WORKDIR
# parallel R: submit job with one MPI master
#mpirun -np 1 R --slave < 00_preparation.R
mpirun -np 1 R --slave < 01_synth_precip.R
mpirun -np 1 R --slave < 02_synth_mle_stat.R
mpirun -np 1 R --slave < 03_synth_mle_nonstat.R
mpirun -np 1 R --slave < 04-spline_synth.R
mpirun -np 1 R --slave < 06-gauge_download.R
mpirun -np 1 R --slave < 07-gauge_mle.R
