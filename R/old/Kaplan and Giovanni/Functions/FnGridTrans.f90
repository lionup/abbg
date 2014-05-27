REAL(8) FUNCTION  FnGridTrans(lx)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lx
REAL(8)             :: lwidth,legrid(ngpe),ledist(ngpe),lvar,ltemp1,ltemp2,ltemp3,ltemp4
INTEGER             :: ij,it


! get boundaries and fill in with equally spaced points
legrid(1)    = -lx*sqrt(Veps)
legrid(ngpe) = lx*sqrt(Veps)
lwidth = (legrid(ngpe)-legrid(1))/real(ngpe-1)
DO ij = 2, ngpe-1
    legrid(ij) = legrid(1) + lwidth*(ij-1)
END DO

! fill in probabilities using normal distribution
CALL cumnor( (legrid(1)+0.5*lwidth)/sqrt(Veps), ledist(1), ltemp2)	
DO ij = 2,ngpe-1
    CALL cumnor ((legrid(ij)+0.5*lwidth)/sqrt(Veps), ltemp1, ltemp2)
    CALL cumnor ((legrid(ij)-0.5*lwidth)/sqrt(Veps), ltemp3, ltemp4)
    ledist(ij) = ltemp1 - ltemp3
END DO
CALL cumnor ((legrid(ngpe)-0.5*lwidth)/sqrt(Veps), ltemp3, ledist(ngpe))
ledist = ledist/sum(ledist)

!find variance
lvar = DOT_PRODUCT(legrid**2,ledist) - DOT_PRODUCT(legrid,ledist)**2

!moment
FnGridTrans  = (Veps -lvar)**2

!put values in globals
DO it = 1,Twork
    edist(it,:) = ledist
    egrid(it,:) = legrid
END DO

END FUNCTION FnGridTrans