rm(list=ls())
#setwd('~/git/abbg2/R/')
source('run.model_nl_st.r')

res <- runMPI(p)
cat('\n',res,'\n')
source('fun.sim.data.r')
sim <- sim.small.sample(p)
save(sim,file='sim_200_399.dat')

