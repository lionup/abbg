SUBROUTINE FnGrossInc(lx,lnet,lf,ldf)
!lx is gross, lnet is net

USE Parameters
USE Globals

IMPLICIT NONE

REAL(8), INTENT(IN) :: lx,lnet
REAL(8), INTENT(OUT) :: lf,ldf


lf = (1.0-pentax - btax)*lx + btax*(lx**(-ptax) + stax)**(-1.0/ptax) - lnet
ldf = (1.0-pentax - btax) + btax* (lx**(-ptax-1.0)) * ((lx**(-ptax) + stax)**(-1.0/ptax - 1.0)) 

END SUBROUTINE FnGrossInc

