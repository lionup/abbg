rm(list = ls())

require(Hmisc)
require(data.table)
require(ggplot2)

source('~/git/abbg/R2/fun.model.solver.nl.ir.r')
source('~/git/abbg/R2/inc.modelsolver.nl.mpi.r')
set.seed(123)


# SETTIG PARAMETERS
p <- list()

# GRIDS DIMENSION - STATE VARIABLES
p$ngpe = 7 #19 #7 			    #transitory component
p$ngpz = 100 #100 #40 #11 			    #permanent component
p$ngpa = 50 #100#50 		    #asset points
p$ngpm = 19  			    #average earnings points
#p$ngpp = p$ngpm * p$ngpz * p$ngpe      #pension points

#DEMOGRAPHIC PARAMETERS PARAMETERS
p$Twork    = 35		#working years 25~59
p$Tret     = 35        #retirement years #60~94
p$Ttot     = p$Twork + p$Tret        #total years

#IMPULSE RESPONSE
p$irb       = 11  #before shock 25~35 years old

#TARGET MOMENTS
p$targetKY = 2.5
p$targetTaxToLabinc = 0.25 #0.17 
p$targetSSAvReplacement = 0.45 

#PARAMETERS FOR GRID CONSTRUCTION
p$pexpgrid = 0.18       #approaches linear as goes to 0, approaches L shaped as goes to Inf
p$amax  = 300000

#SIMULATION PARAMETERS	
p$nsim = 50000 #5000
p$N = 999999

#EARNINGS PROCESS
p$Veps =  0.05
p$Vz0  =  0.15
p$rho  =  1
p$delta =  0.2
p$Veta_rho1 =  0.01   #Veta if rho==1
p$tau   = 0.15

#INTEREST RATE
p$R = 1.03      #annual gross interest rate

#PREFERENCE PARMETERS
p$gam =   2
p$bet =   0.9678168 # 0.961165 #zbl #1/p$R

#BORROWING LIMIT: SET TO VERY LARGE NEGATIVE NO. FOR NBL
p$borrowlim = -100000000.0 #0.0 #-100000000.0

#GOVERNMENT PARAMETERS
#gouveia strauss 
#p$stax = 2.0e-4   #guess: is chosen optimally
p$ptax = 0.768
p$btax = 0.258

#other taxes and benefits
p$cfloor = -100000000000.0  
p$pentax  = 0.0    #payroll tax   
p$rtax = 0.0 #tax on interest income
p$otax = 0.0 #tax on old age pensions, check option in Parameters
p$pencapfrac = 2.2 #cap on (pre-tax) earnings that contribute to pension index, 
p$Rnet = 1.0 + (1.0-p$rtax)*(p$R-1)      #annual after tax interest rate

#OPTIONS
p$Display  = 1
p$mode <- 'multicore' #'serial' #'multicore' #'mpi' 

start_time = proc.time()[3]  
models  <- comp.solveModel(p)
save( models, file = paste('sim.ir.rule.nbl',p$ngpz,'dat',sep='.') ) 
cat(paste('\ntime for decision rules: ' , proc.time()[3] -  start_time ))

#load('~/git/abbg/R2/sim.ir.rule.100.dat')

triquant <- c(0.1,0.5,0.9)

for( tau0 in triquant ){
	samples  <- comp.samples(models, tau0)
	save( samples, file=paste('sim.ir.sample.nbl',p$ngpz,tau0,'dat',sep='.') ) 
	#load( paste('~/git/abbg/R2/sim.ir.sample',p$ngpz,tau0,'dat',sep='.') )

	for( tau1 in triquant ){
		start_time = proc.time()[3]  
		moments  <- comp.moments(models, samples, tau0, tau1)
		save( moments, file=paste('sim.ir.res.nbl',p$ngpz,tau0,tau1,'dat',sep='.') ) 
		cat(paste('\ntime for tau0/tau1:', tau0, '/', tau1, 'is', proc.time()[3] -  start_time ))
	}
}


