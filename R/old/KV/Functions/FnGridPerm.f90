REAL(8) FUNCTION  FnGridPerm(lx)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN)         :: lx
REAL(8),DIMENSION(Twork)    :: lwidth,lvar
REAL(8)                     :: ltemp1,ltemp2,ltemp3,ltemp4
INTEGER             :: iz1,iz2,it


! get boundaries and fill in with equally spaced points
DO it = 1,Twork
    zgrid(it,1)    = -lx*sqrt(varz(it))
    zgrid(it,ngpz) = lx*sqrt(varz(it))
    lwidth(it) = (zgrid(it,ngpz)-zgrid(it,1))/real(ngpz-1)
    DO iz1 = 2, ngpz-1
        zgrid(it,iz1) = zgrid(it,1) + lwidth(it)*(iz1-1)
    END DO        
END DO

! fill in transition matrix using normal distribution
DO it = 1,Twork-1
    DO iz1 = 1,ngpz
        CALL cumnor ((zgrid(it+1,1)+0.5*lwidth(it+1)-rho*zgrid(it,iz1))/sqrt(Vetavec(it)), ztrans(it,iz1,1), ltemp2)
	    DO iz2 = 2,ngpz-1
	        CALL cumnor ((zgrid(it+1,iz2)+0.5*lwidth(it+1)-rho*zgrid(it,iz1))/sqrt(Vetavec(it)), ltemp1, ltemp2)
	        CALL cumnor ((zgrid(it+1,iz2)-0.5*lwidth(it+1)-rho*zgrid(it,iz1))/sqrt(Vetavec(it)), ltemp3, ltemp4)
            ztrans(it,iz1,iz2) = ltemp1 - ltemp3
        END DO
        CALL cumnor ((zgrid(it+1,ngpz)-0.5*lwidth(it+1)-rho*zgrid(it,iz1))/sqrt(Vetavec(it)), ltemp3,ztrans(it,iz1,ngpz))
        ztrans(it,iz1,:) = ztrans(it,iz1,:)/sum(ztrans(it,iz1,:))
    END DO
END DO

!find distribution at first period
CALL cumnor( (zgrid(1,1)+0.5*lwidth(1))/sqrt(Vz0), zdist(1,1), ltemp2)	
DO iz1 = 2,ngpz-1
    CALL cumnor ((zgrid(1,iz1)+0.5*lwidth(1))/sqrt(Vz0), ltemp1, ltemp2)
    CALL cumnor ((zgrid(1,iz1)-0.5*lwidth(1))/sqrt(Vz0), ltemp3, ltemp4)
    zdist(1,iz1) = ltemp1 - ltemp3
END DO
CALL cumnor ((zgrid(1,ngpz)-0.5*lwidth(1))/sqrt(Vz0), ltemp3, zdist(1,ngpz))
zdist(1,:) = zdist(1,:)/sum(zdist(1,:))

!find unconditional distributions
DO it = 2,Twork
    zdist(it,:) = MATMUL(zdist(it-1,:),ztrans(it-1,:,:))
    zdist(it,:) = zdist(it,:) /sum(zdist(it,:) )
END DO


!find variance
DO it = 1,Twork
    lvar(it) = DOT_PRODUCT(zgrid(it,:)**2,zdist(it,:)) - DOT_PRODUCT(zgrid(it,:),zdist(it,:))**2
END DO

!moment
FnGridPerm  = SUM((varz -lvar)**2)

!put values in globals
varzapprox = lvar

END FUNCTION FnGridPerm