#!/bin/bash

echo "startin qsub script file"
source ~/.bash_profile


# here's the SGE directives
# ------------------------------------------
#$ -q batch.q   		 # <- the name of the Q you want to submit to
#$ -pe mpich 51 		 #  <- load the openmpi parallel env w/ 4 slots
#$ -S /bin/bash  		 # <- run the job under bash
#$ -N abbg      		 # <- name of the job in the qstat output
#$ -o mpi_abbg.out   # <- name of the output file.
#$ -e mpi_abbg.out   # <- name of the stderr file.
#$ -cwd

module load openmpi
module load open64
module load r/3.0.1

echo "calling mpirun now"

mpirun -np 51 /data/uctprgu/R/x86_64-unknown-linux-gnu-library/3.0/snow/RMPISNOW -q < run.model_nl_kg.r > mpi_abbg.rout 

## call via: qsub qsub_start_mpi.sh 

