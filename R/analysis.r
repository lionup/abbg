
ssh uctprgu@hpc.econ.ucl.ac.uk

scp uctprgu@hpc.econ.ucl.ac.uk:/data/uctprgu/git/abbg2/R/cohort* .

scp US.dat uctprgu@hpc.econ.ucl.ac.uk:/data/uctprgu/git/lg2/Stata/result3/cmes2014/


rm(list = ls())
setwd('~/git/lgbx/R')
load('evaluations.dat')
I = which(param_data$value == min(param_data$value,na.rm=TRUE))[[1]]
pp = c(param_data[I,])

require(plyr)
list2df <- function(ll) {
  return(ldply(ll,function(l){ return(data.frame(rbind(unlist(l))))}))
}

sm <- list2df(pp)
names(sm) <- c('names', 'value')
sm$value <- round(sm$value,4)
sm

params_to_sample  = cf$params_to_sample
params_to_sample2 = paste('p',cf$params_to_sample,sep='.')
require(data.table)
param_data = data.table(param_data) 
p = cf$initial_value
p[params_to_sample] = param_data[I,params_to_sample2,with=FALSE]



load('maxres.dat')

load(mcf$file_chain)
p.mean=lapply(param_data[,mcf$params_to_sample],mean,na.rm=TRUE)
p.mean=c(p.mean,mcf$initial_value[setdiff(names(mcf$initial_value),names(p.mean))])
meanres=MOPT_OBJ_FUNC(p.mean)


source('~/git/Utils/R/inc.utils.r')



load('param_error.dat')
p<-as.list( head(tail(per,2),1) )
save(p,file='p.dat')
