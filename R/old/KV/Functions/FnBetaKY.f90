REAL(8) FUNCTION FnBetaKY(lbet)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lbet
REAL(8)             :: lcapital,lincome

bet = lbet
write(*,*) ' Beta guess: ',bet

CALL Decisions
CALL Simulate

lcapital = SUM(SUM(asim(:,1:Ttot), DIM=1)*popsize)
IF (UseNetIncKY==1) THEN
    lincome = (R-1)*lcapital + SUM(SUM(ysim(:,:), DIM=1)*popsize) 
ELSE
    lincome = (R-1)*lcapital + SUM(SUM(ypresim(:,:), DIM=1)*popsize) 
END IF

simKY =lcapital/lincome
FnBetaKY = simKY-targetKY

write(*,*) ' Simulated - Target KY: ',FnBetaKY


END FUNCTION FnBetaKY
