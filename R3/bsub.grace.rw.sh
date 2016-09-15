#!/bin/bash -l

#BSUB -q cowles
##BSUB -q shared
#BSUB -W 1:00
#BSUB -J abbg_R3

#BSUB -n 20
#BSUB -R "span[hosts=1]"
##BSUB -n 40
##BSUB -R "span[ptile=20]"

#BSUB -o grace.sys.%J   # output file name in which %J is replaced by the job ID
#BSUB -N                   # sent email

# 9. load modules
echo "Slot count: $LSB_DJOB_NUMPROC"
echo "Slot distribution: $LSB_MCPU_HOSTS"

module load Apps/R/3.2.2-generic
module load Libs/GSL/1.16

R --no-save --slave < ~/git/abbg/R3/run.model.rw.mpi.r > grace.$LSB_JOBID
