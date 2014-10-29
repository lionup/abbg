cat(' Beta guess: ',lbet, '\n')

moments  <- comp.solveModel(p)

lcapital = SUM(SUM(asim[,1:Ttot], DIM=1)*popsize)
lincome = (R-1)*lcapital + SUM(SUM(ypresim, DIM=1)*popsize) 

simKY =lcapital/lincome
FnBetaKY = simKY-targetKY

cat(' Simulated - Target KY: ',FnBetaKY, '\n')
