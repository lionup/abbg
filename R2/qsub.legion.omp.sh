#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=1:0:0
#$ -l mem=8G
#$ -N abbg
#$ -l thr=12
#$ -wd /home/uctprgu/Scratch/abbg
#$ -j y
#$ -M uctprgu@ucl.ac.uk	# send notifications to
#$ -m e # notify about end of job

# 9. load modules
module unload compilers
module load recommended/r-new
#module load gsl/1.15/gnu.4.6.3

#R --no-save --slave < $HOME/git/abbg/R2/run.model.rw.mpi.r > r.output.$JOB_ID
R --no-save --slave < $HOME/git/abbg/R2/run.model.nl.mpi.r > r.output.$JOB_ID
