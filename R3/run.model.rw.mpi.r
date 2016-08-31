rm(list = ls())
require(Hmisc)
require(data.table)
require(ggplot2)
require(EQL)

setwd('~/git/abbg/R3')
source('~/git/abbg/R3/fun.model.solver.rw.mpi.r')
set.seed(123)

# SETTIG PARAMETERS
p <- list()

# GRIDS DIMENSION - STATE VARIABLES
p$ngpe = 7 #19 #7 			    #transitory component
p$ngpz = 39 #11 			    #permanent component
p$ngpa = 50#100#50 		    #asset points
p$ngpm = 19  			    #average earnings points
#p$ngpp = p$ngpm * p$ngpz * p$ngpe      #pension points

#DEMOGRAPHIC PARAMETERS PARAMETERS
p$age_min = 25
p$age_re  = 61
p$age_max = 93
p$Twork  = (p$age_re - p$age_min)/2
p$Tret   = (p$age_max - p$age_re)/2 +1 #periods after retirement
p$Ttot   = p$Twork + p$Tret        #total years

#TARGET MOMENTS
p$targetKY = 2.5
p$targetTaxToLabinc = 0.25 #0.17
p$targetSSAvReplacement = 0.45

#PARAMETERS FOR GRID CONSTRUCTION
p$pexpgrid = 0.18       #approaches linear as goes to 0, approaches L shaped as goes to Inf
p$amax  = 300000

#SIMULATION PARAMETERS
p$nsim = 50000 #5000
#p$N = 999999

#EARNINGS PROCESS
p$Veps =  0.0887
p$Vz0  =  0.1127
p$rho  =  1
#p$delta =  0.2
p$Veta_rho1 =  0.0183   #Veta if rho==1
#p$tau   = 0.15

#INTEREST RATE
p$R = 1.03^2      #annual gross interest rate

#PREFERENCE PARMETERS
p$gam =   2 #2 #15
p$bet =   1/p$R

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
moments  <- comp.solveModel(p)
cat(paste('\ntotal seconds to solve the program: ' , proc.time()[3] -  start_time ))

save(p, moments, file='rw_nbl_full.dat')
