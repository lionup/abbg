SUBROUTINE Simulate

USE Parameters
USE Globals
USE Procedures
USE random

IMPLICIT NONE

INTEGER		:: iseed(1)

INTEGER		:: in,it,itemp,lim(2)
REAL(8)	    :: ltemp,lue,lm,lmu,lsig2,lp(2)


!set random number seed
iseed(1) = 6969
CALL RANDOM_SEED(put = iseed)   

!initial earnings
CALL RandomDiscrete(nsim, zsimI(:,1), ngpz, zdist(1,:))
CALL RandomDiscrete(nsim, esimI(:,1), ngpe, edist(1,:))
zsim(:,1) = zgrid(1,zsimI(:,1))
esim(:,1) = egrid(1,esimI(:,1))
DO in = 1,nsim
    ypresim(in,1) = ypregrid(1,zsimI(in,1),esimI(in,1))
    yavsim(in,1) = min(ypresim(in,1),pencap)
    CALL FindLinProb1(yavsim(in,1),mgrid(1,:),lim,lp)
    CALL RandomDiscrete1(itemp,2,lp)
    msimI(in,1) = lim(itemp)
    msim(in,1) = mgrid(1,msimI(in,1))
    ysim(in,1) = ygrid(1,zsimI(in,1),esimI(in,1))
END DO

!initial assets
IF (IntitialWealthDist==1) THEN
    DO in = 1,nsim
        CALL RandomDiscrete1(itemp, 75, initwealthdist(2,:))
        asim(in,1) = initwealthdist(1,itemp)*ysim(in,1)
        asim(in,1) = max(asim(in,1), agrid(1,1))
    END DO
ELSE
    asim(:,1) = 0.0
END IF

!do t=1 separately
it = 1
DO in = 1,nsim
    trsim(in,it)    = max( cfloor - (Rnet*asim(in,it)+ ysim(in,it) - agrid(it+1,1)) ,0.0)
    CALL LinInterp1 (ngpa,agrid(it,:),con(it,:,msimI(in,it),zsimI(in,it),esimI(in,it)),asim(in,it)+trsim(in,it),csim(in,it))
    IF (trsim(in,it)>0.0) csim(in,it) = min(cfloor, csim(in,it))
    xsim(in,it)     = Rnet*asim(in,it)+ ysim(in,it) + trsim(in,it)
    asim(in,it+1)   = xsim(in,it) - csim(in,it)
    IF (asim(in,it+1)<agrid(it+1,1))THEN
        asim(in,it+1) = agrid(it+1,1)
        csim(in,it) = xsim(in,it) - asim(in,it+1)
    END IF        
    IF(isnan(csim(in,it))) THEN
        write(*,*) 'bad consumption'
    END IF
    
END DO

!working life
DO it = 2,Twork
    DO in = 1,nsim
        
        CALL RandomDiscrete1(zsimI(in,it), ngpz, ztrans(it-1,zsimI(in,it-1),:))
        CALL RandomDiscrete1(esimI(in,it), ngpe, edist(it,:))
        
        zsim(in,it) = zgrid(it,zsimI(in,it))
        esim(in,it) = egrid(it,esimI(in,it))
        ysim(in,it) = ygrid(it,zsimI(in,it),esimI(in,it))
        ypresim(in,it) = ypregrid(it,zsimI(in,it),esimI(in,it))
        yavsim(in,it) = ((it-1)*msim(in,it-1) + min(ypresim(in,it),pencap))/real(it)
        CALL FindLinProb1(yavsim(in,it),mgrid(it,:),lim,lp)
        CALL RandomDiscrete1(itemp,2,lp)
        msimI(in,it) = lim(itemp)
        msim(in,it) = mgrid(it,msimI(in,it))
        
        IF (it<Twork) THEN
            trsim(in,it)    = max( cfloor - (Rnet*asim(in,it)+ ysim(in,it) - agrid(it+1,1)) ,0.0)
            CALL LinInterp1 (ngpa,agrid(it,:),con(it,:,msimI(in,it),zsimI(in,it),esimI(in,it)),asim(in,it),csim(in,it))
        ELSEIF (it==Twork) THEN
            psimI(in) = pind(msimI(in,Twork),zsimI(in,Twork),esimI(in,Twork))
            trsim(in,it)    = max( cfloor - (Rnet*asim(in,it)+ ysim(in,it) - agridret(1,1,psimI(in)) ) ,0.0)
            CALL LinInterp1 (ngpa,agrid(it,:),con(it,:,msimI(in,it),zsimI(in,it),esimI(in,it)),asim(in,it),csim(in,it))
        END IF
        xsim(in,it)     = Rnet*asim(in,it)+ ysim(in,it) + trsim(in,it)
        asim(in,it+1)   = xsim(in,it) - csim(in,it)
        
        IF (it<Twork) THEN    
            IF (asim(in,it+1)<agrid(it+1,1))THEN
                asim(in,it+1) = agrid(it+1,1)
                csim(in,it) = xsim(in,it) - asim(in,it+1)
            END IF        
        ELSEIF (it==Twork) THEN
            IF (asim(in,it+1)<agridret(1,1,psimI(in)))THEN
                asim(in,it+1) = agridret(1,1,psimI(in))
                csim(in,it) = xsim(in,it) - asim(in,it+1)
            END IF                
        END IF
        
        IF(isnan(csim(in,it))) THEN
            write(*,*) 'bad consumption'
        END IF
    
    END DO
END DO

!retirement
DO it = 1,Tret
    DO in = 1,nsim
        ysim(in,Twork+it) = pgrid(it,psimI(in))
        ypresim(in,Twork+it) = ppregrid(it,psimI(in))
        IF (it<Tret) trsim(in,Twork +it)    = max( cfloor - (Rnet*asim(in,it)+ ysim(in,it) - agridret(it+1,1,psimI(in))) ,0.0)        
        CALL LinInterp1 (ngpa,agridret(it,:,psimI(in)),conret(it,:,psimI(in)),asim(in,Twork+it),csim(in,Twork+it))
        xsim(in,Twork+it)     = Rnet*asim(in,Twork+it)+ ysim(in,Twork+it)  + trsim(in,Twork +it)
        asim(in,Twork+it+1) = (xsim(in,Twork+it)  - csim(in,Twork+it))/annprem(it)
        
        IF (it<Tret) THEN
            IF (asim(in,Twork+it+1)<agridret(it+1,1,psimI(in)))THEN
                asim(in,Twork+it+1) = agridret(it+1,1,psimI(in))
                csim(in,Twork+it) = xsim(in,Twork+it) - asim(in,Twork+it+1)*annprem(it)
            END IF        
        ELSEIF (it==Tret) THEN
            IF (asim(in,Twork+it+1) .ne. 0)THEN
                asim(in,Twork+it+1) = 0
                csim(in,Twork+it) = xsim(in,Twork+it)
            END IF             
        END IF
        
    END DO
END DO

END SUBROUTINE Simulate