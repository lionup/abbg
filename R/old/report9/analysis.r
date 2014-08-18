
ssh uctprgu@hpc.econ.ucl.ac.uk

scp uctprgu@hpc.econ.ucl.ac.uk:/data/uctprgu/git/abbg/R/cohort30.dat .

scp US.dat uctprgu@hpc.econ.ucl.ac.uk:/data/uctprgu/git/lg2/Stata/result3/cmes2014/

source('fun.sim.data.r')
load('sim_100_99.dat')
#simdata <-data.matrix(sim)
#save(simdata,file='simdata.dat')
#require(R.matlab)
#writeMat('simdata.mat',simdata=simdata)

#sim$t <-1:p$T
#sim_y <- sim[,c("id","t","Y"),with=F]
#wide_y <- reshape(sim_y, idvar='id', timevar='t', direction='wide') #each person has all age in a row
#Y <-data.matrix(wide_y[,-1,with=F])
#writeMat('Y.mat',Y=Y)
persis <- sim.persis(p,sim)
summary(c(persis))

