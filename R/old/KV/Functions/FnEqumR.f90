REAL(8) FUNCTION FnEqumR(lR)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lR

R = lR
Rnet = 1.0 + (1.0-rtax)*(R-1)   
write(*,*) '  R guess: ',R

CALL Grids
CALL Decisions
CALL Simulate


FnEqumR = SUM(SUM(asim(:,1:Ttot), DIM = 1)*popsize)/SUM(SUM(ysim(:,:), DIM = 1)*popsize)

write(*,*) '  Net assets relative to income: ',FnEqumR


END FUNCTION FnEqumR
