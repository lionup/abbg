#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=3:0:0
#$ -l mem=8G
#$ -N abbg-ser
#$ -wd /home/uctprgu/Scratch/abbg
#$ -j y
#$ -M uctprgu@ucl.ac.uk	# send notifications to
#$ -m e # notify about end of job

module unload compilers
module load recommended/r-new

R --no-save --slave < $HOME/git/abbg/R2/run.model.nl.ir.r >r.output.$JOB_ID

