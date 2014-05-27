SUBROUTINE LoadData

USE Parameters
USE Globals

IMPLICIT NONE

INTEGER :: it
REAL(8) :: unconsurprob(Tret)

!log experience profile
OPEN(1, FILE = InputDir // 'kappasmooth.txt')
READ(1,*) kappa
CLOSE(1)
kappa = kappa + log(0.75)

!survival probabilities
OPEN(1, FILE = InputDir // 'surprobsmooth.txt')
READ(1,*) surprob
CLOSE(1)

!age-specific variances
IF(AgeSpecificVariances==1) THEN
	OPEN(1, FILE = InputDir // 'varetaage.txt')
	READ(1,*) varetaage
	CLOSE(1)
	Vetavec = Veta*varetaage
END IF

unconsurprob(1) = 1.0 
DO it = 1,Tret
    IF(it>1) unconsurprob(it) = unconsurprob(it-1)* surprob(it-1)
	IF (AnnuityMarkets == 1) THEN
		annprem(it) = surprob(it)
	ELSE
		annprem(it) = 1.0
	END IF
END DO

DO it = 1,Twork
    popsize(it) = 1.0/real(Twork + sum(unconsurprob))
END DO
DO it = 1,Tret
    popsize(Twork+it) = unconsurprob(it)/real(Twork + sum(unconsurprob))
END DO
    
!initial wealth distribution
OPEN(1, FILE = InputDir // 'initwealthdist.txt')
READ(1,*) initwealthdist
CLOSE(1)
    
END SUBROUTINE LoadData