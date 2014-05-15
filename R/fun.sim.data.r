##data for stephane
sim.origin.sample <- function(p){
  load('even.dat')
  age_full = seq(p$age_min, p$age_max, 2)
  last_age = p$age_re - 2*6  #oldest cohort ini age
  age_nore = seq(p$age_min, last_age, 2)
  attach(moments)
  income <- t(ytList)
  consumption <- t(ctList)
  saving     <- t(stList)
  eta <- t(etaList)
  eps <- t(epsList)
  detach(moments)
  require(data.table)
  sim_data <- data.table(id=1:p$nsim, age=rep(age_full,each=p$nsim),
              income=c(income),consumption=c(consumption),
              saving=c(saving),eta=c(eta),eps=c(eps))
  setkey(sim_data,id)

  #for each age, randomly pick 40 people
  sub_id <- array( sample(1:p$nsim, length(age_nore)*40), 
    dim=c(length(age_nore),40) )

  sim <- data.table()
  for( t in 1:length(age_nore) ){
    sim_sub <- sim_data[J(sub_id[t,])]
    sim_sub <- subset(sim_sub, age>=age_nore[t] & age<=age_nore[t]+10)
    sim <- rbind(sim,sim_sub)
  }
  #######################################
  p$age_min = 31
  p$age_re  = p$age_min+36 #first period income drop almost half
  p$age_max = p$age_min+50

  load('odd.dat')
  age_full = seq(p$age_min, p$age_max, 2)
  last_age = p$age_re - 2*6
  age_nore = seq(p$age_min, last_age, 2)
  attach(moments)
  income <- t(ytList)
  consumption <- t(ctList)
  saving     <- t(stList)
  eta <- t(etaList)
  eps <- t(epsList)
  detach(moments)
  sim_data2 <- data.table(id=(p$nsim+1):(2*p$nsim), age=rep(age_full,each=p$nsim),
              income=c(income),consumption=c(consumption),
              saving=c(saving),eta=c(eta),eps=c(eps))
  setkey(sim_data2,id)

  #for each age, randomly pick 40 people
  sub_id <- array( sample((p$nsim+1):(2*p$nsim), length(age_nore)*40), 
    dim=c(length(age_nore),40) )

  for( t in 1:length(age_nore) ){
    sim_sub <- sim_data2[J(sub_id[t,])]
    sim_sub <- subset(sim_sub, age>=age_nore[t] & age<=age_nore[t]+10)
    sim <- rbind(sim,sim_sub)
  }
  setkey(sim,id,age)

  simdata <-data.matrix(sim)

  p$age_min = 30
  p$age_re  = p$age_min+36 #first period income drop almost half
  p$age_max = p$age_min+50
  
  return(simdata)
}