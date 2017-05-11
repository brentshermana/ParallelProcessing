#!/bin/bash
#SBATCH --job-name TestJob
##SBATCH --nodes 1
#SBATCH --ntasks 16
##SBATCH --mem 100MB
##SBATCH --mem-per-cpu 10M
#SBATCH --partition class
#SBATCH --time 00:05:00
#SBATCH --error=TestJob.%J.err
#SBATCH	--output=TestJob.%j.out
######################################
######## END OF SLURM OPTIONS ########
######################################

## From this point on, simply specify what you want
## your job to do.

module load mpich
module load gcc

mpiexec ~/BrentStuff/MPI/cpi2

