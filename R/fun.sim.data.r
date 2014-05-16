##data for stephane
sim.origin.sample <- function(p){
  load('even.dat')
  age_full = seq(p$age_min, p$age_max, 2)
  last_age = p$age_re - 2*6  #oldest cohort ini age
  age_nore = seq(p$age_min, last_age, 2)

  require(data.table)
  sim_data <- with( moments, data.table(id=rep(1:p$nsim,each=p$nage+p$ntr), 
    age=age_full, income=c(ytList),consumption=c(ctList),
    saving=c(stList),eta=c(etaList),eps=c(epsList)) )
  sim_data[,Y:=eta+eps]
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
  load('odd.dat')
  age_full = seq(p$age_min, p$age_max, 2)
  last_age = p$age_re - 2*6
  age_nore = seq(p$age_min, last_age, 2)

  sim_data2 <- with( moments, data.table(id=rep((p$nsim+1):(2*p$nsim),each=p$nage+p$ntr), 
    age=age_full, income=c(ytList),consumption=c(ctList),
    saving=c(stList),eta=c(etaList),eps=c(epsList)) )
  sim_data2[,Y:=eta+eps]
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
  return(sim)
}

