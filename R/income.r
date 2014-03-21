require(R.matlab)
setwd('~/git/abbg/R')
data <- readMat('data_hermite.mat')

# two periods
#declare parameters
p <- list()
p$nw <- 11  
p$nt <- 31  #from 30 to 60

xw <- array(0, dim=c(p$nt, p$nw))
y1 <- dataResqfinal.e0[1,]

#get the quantile
Vectau <- 1/(1+p$nw)

#	create the grid nodes of eta_1
# Resqtrue_e0 is the 
MatAGE1 <- c(1,0)


#	create the max and min nodes on the grid of eta_2
# then discretize the grid into ne points

#for each grid node of eta_1, locate the quantiles for all grid nodes of eta_2

#utting the unit interval into 11 portions by the quantiles calculated 

#the transition probability is the length of each protions
