rm(list=ls())
setwd('~/git/abbg2/R/')
source('run.model_nl_pa.r')
source('fun.model_solver_nl_st.r')
source('inc.modelsolver_nl_st.r')

res = runMPI(p)

