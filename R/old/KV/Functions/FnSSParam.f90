REAL(8) FUNCTION FnSSParam(lsspar)
USE Parameters
USE Globals
USE Procedures
IMPLICIT NONE

REAL(8), INTENT(IN) :: lsspar
INTEGER         :: it,iz,ie,ip,im,in,lmloc(2)
REAL(8)     ::ltotpen,lmearns,lpsim(nsim),lp(2),lavearns,lmpen
REAL(8), EXTERNAL        :: FnTax

IF (Display==1)write(*,*) ' SS benefit parameter: ', lsspar

IF (BendPointsPostTax==1) THEN
    lavearns = SUM(avearnspost)/real(Twork) 
ELSE
    lavearns = SUM(avearnspre)/real(Twork) 
END IF

IF (ScaleBendPoints==1) THEN    
    ssbend1 = 0.18*lsspar*lavearns 
    ssbend2 = 1.1*lsspar*lavearns 
ELSE
    ssbend1 = 0.18*lavearns
    ssbend2 = 1.1*lavearns
END IF
 
ltotpen = 0.0
ip = 1
DO im = 1,ngpm
DO iz = 1,ngpz
DO ie = 1,ngpe
    pind(im,iz,ie) = ip    
    DO it = 1,Tret
        IF (UseFinalZPension==1 )THEN
            lmearns = dexp(kappa(Twork) + zgrid(Twork,iz))
        ELSE
            !lmearns = ((Twork-1)*mgrid(Twork-1,im) + min(ypregrid(Twork,iz,ie),pencap))/real(Twork)
            lmearns = mgrid(Twork,im)
        END IF
        
        IF (lmearns<= ssbend1) THEN
            ppregrid(it,ip) = 0.9*lmearns
        ELSEIF (lmearns<= ssbend2) THEN
            ppregrid(it,ip) = 0.9*ssbend1 + 0.32*(lmearns - ssbend1)
        ELSE
            ppregrid(it,ip) = 0.9*ssbend1 + 0.32*(ssbend2-ssbend1) + 0.15*(lmearns - ssbend2)
        END IF

        IF (ScaleBendPoints==0) ppregrid(it,ip) = ppregrid(it,ip)*lsspar
        !pgrid(it,ip) = ppregrid(it,ip)*(1.0-otax)
        pgrid(it,ip) = ppregrid(it,ip) - FnTax(ppregrid(it,ip)*0.85)
    
    END DO        
    ip = ip+1
END DO
END DO
END DO


IF (MatchAggPensionBen==1) THEN

    IF (UseFinalZPension==1) THEN

        ltotpen = 0.0
        DO iz = 1,ngpz
        DO ie = 1,ngpe
            ip = pind(1,iz,ie)
            DO it = 1,Tret
                ltotpen = ltotpen + ppregrid(it,ip)*zdist(Twork,iz)*edist(Twork,ie)*popsize(Twork+it)
            END DO        
            ip = ip+1
        END DO
        END DO


    ELSE
        DO in = 1,nsim
            CALL FindLinProb1(yavsim(in,Twork-1),mgrid(Twork-1,:),lmloc,lp)
            CALL RandomDiscrete1(im,2,lp)
            lpsim(in) = ppregrid(1,pind(lmloc(im),zsimI(in,Twork),esimI(in,Twork)))
        END DO

        ltotpen = 0.0
        DO it = 1,Tret
            ltotpen = ltotpen + popsize(Twork+it)*sum(lpsim)/real(nsim)
        END DO
END IF

    IF (UseNetIncToScaleSS == 1) THEN
        FnSSParam = ltotpen/totlabincpost - targetSSBenToLabinc
        IF (Display==1)write(*,*) ' SS Ben / Post-tax labor income: ', (ltotpen/totlabincpost)*100, '%'
    ELSE
        FnSSParam = ltotpen/totlabincpre - targetSSBenToLabinc
        IF (Display==1)write(*,*) ' SS Ben / Pre-tax labor income: ', (ltotpen/totlabincpre)*100, '%'
    END IF
ELSE !match benefits based on guy with average AIME
    lmearns =  sum(min(avearnspre(:),pencap))/real(Twork)
    IF (lmearns<= ssbend1) THEN
        lmpen = 0.9*lmearns
    ELSEIF (lmearns<= ssbend2) THEN
        lmpen = 0.9*ssbend1 + 0.32*(lmearns - ssbend1)
    ELSE
        lmpen = 0.9*ssbend1 + 0.32*(ssbend2-ssbend1) + 0.15*(lmearns - ssbend2)
    END IF
    IF (ScaleBendPoints==0) lmpen = lmpen*lsspar
    
    FnSSParam = lmpen/(SUM(avearnspre)/real(Twork)) - targetSSAvReplacement
    IF (Display==1)write(*,*) ' SS Ben for Mean AIME/ Av Pre-tax Earns: ', (lmpen/(SUM(avearnspre)/real(Twork)))*100, '%'

END IF
    
END FUNCTION FnSSParam

