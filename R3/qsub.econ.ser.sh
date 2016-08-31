#!/bin/bash

echo "startin qsub script file"
source ~/.bash_profile
date

# here's the SGE directives
# ------------------------------------------
#$ -q batch.q   		 # <- the name of the Q you want to submit to
#$ -S /bin/bash  		 # <- run the job under bash
#$ -N js_ser         		 # <- name of the job in the qstat output
#$ -o js_ser.out   	     # <- name of the output file.
#$ -j y
#$ -cwd

echo "loaded modules"
module load sge/2011.11
module add r/3.1.3
module list

echo "calling R now"
R --no-save -q < run.model.nl.mpi.r > js_ser.Rout
