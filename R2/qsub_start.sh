#!/bin/bash

echo "startin qsub script file"
source ~/.bash_profile


# here's the SGE directives
# ------------------------------------------
#$ -q batch.q   		 # <- the name of the Q you want to submit to
#$ -S /bin/bash  		 # <- run the job under bash
#$ -N abbg      		 # <- name of the job in the qstat output
#$ -o sysinf.$JOB_ID   # <- name of the output file.
#$ -j y
#$ -cwd

module load r

echo "calling mpirun now"

#R --no-save < run.model_abbg.r > r.output.$JOB_ID
R --no-save < run.model_rwR.r > r.output.$JOB_ID

## call via: qsub qsub_start_mpi.sh 

