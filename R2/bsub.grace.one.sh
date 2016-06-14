#!/bin/bash -l

#BSUB -q cowles
#BSUB -W 02:00
#BSUB -J abbg

#BSUB -n 20
#BSUB -R "span[hosts=1]" 
##BSUB -n 40
##BSUB -R "span[ptile=20]" 

#BSUB -o grace.output.%J   # output file name in which %J is replaced by the job ID
##BSUB -N                   # sent email

# 9. load modules
echo "Slot count: $LSB_DJOB_NUMPROC"
echo "Slot distribution: $LSB_MCPU_HOSTS"

module load Apps/R/3.2.2-generic
module load Libs/GSL/1.16

R --no-save --slave < ~/git/abbg/R2/run.cons.reg.st.parallel.r > grace.R.$LSB_JOBID
