#!/bin/bash

echo "startin qsub script file"
source ~/.bash_profile


# here's the SGE directives
# ------------------------------------------
#$ -q batch.q   		 # <- the name of the Q you want to submit to
#$ -pe mpich 40  #  <- load the openmpi parallel env w/ 4 slots
#$ -S /bin/bash  		 # <- run the job under bash
#$ -N abbg      		 # <- name of the job in the qstat output
#$ -o sysinf.$JOB_ID   # <- name of the output file.
#$ -j y
#$ -cwd

module add openmpi/gcc
#module load openmpi
#module load open64
module load gcc
module add r

echo "calling mpirun now"

mpirun -np 40 /data/uctprgu/R/x86_64-unknown-linux-gnu-library/3.1/snow/RMPISNOW -q < run.model_abbg.mpi.r > r.output.$JOB_ID
#R --no-save < run.model_rwR.r > r.output.$JOB_ID

## call via: qsub qsub_start_mpi.sh 

