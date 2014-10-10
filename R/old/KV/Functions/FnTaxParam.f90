REAL(8) FUNCTION FnTaxParam(lstax)
USE Parameters
USE Globals
IMPLICIT NONE

REAL(8), INTENT(IN) :: lstax
INTEGER         :: it,iz,ie
REAL(8)     ::ltotlabincpre,ltotlabincpost,ltottax
REAL(8), EXTERNAL        :: FnTax

stax = lstax
IF (Display==1)write(*,*) ' Tax function parameter: ', stax

ltottax = 0.0
ltotlabincpre = 0.0
ltotlabincpost = 0.0

DO it = 1,Twork
    avearnspre(it) = 0.0
    avearnspre2(it) = 0.0
    avearnspost(it) = 0.0
    avearnspost2(it) = 0.0
    avlearnspre(it) = 0.0
    avlearnspre2(it) = 0.0
    avlearnspost(it) = 0.0
    avlearnspost2(it) = 0.0

    DO iz = 1,ngpz
    DO ie = 1,ngpe
        ypregrid(it,iz,ie) = dexp(kappa(it) +egrid(it,ie) + zgrid(it,iz))
        ygrid(it,iz,ie) = ypregrid(it,iz,ie) - FnTax(ypregrid(it,iz,ie))
        
        ltotlabincpre = ltotlabincpre + ypregrid(it,iz,ie)*zdist(it,iz)*edist(it,ie)*popsize(it)
        ltotlabincpost = ltotlabincpost + ygrid(it,iz,ie)*zdist(it,iz)*edist(it,ie)*popsize(it)
        ltottax = ltottax + (FnTax(ypregrid(it,iz,ie)) - pentax*ypregrid(it,iz,ie))*zdist(it,iz)*edist(it,ie)*popsize(it)
        avearnspre(it) = avearnspre(it) + ypregrid(it,iz,ie)*zdist(it,iz)*edist(it,ie)
        avearnspre2(it) = avearnspre2(it) + (ypregrid(it,iz,ie)**2)*zdist(it,iz)*edist(it,ie)
        avearnspost(it) = avearnspost(it) + ygrid(it,iz,ie)*zdist(it,iz)*edist(it,ie)
        avearnspost2(it) = avearnspost2(it) + (ygrid(it,iz,ie)**2)*zdist(it,iz)*edist(it,ie)
        avlearnspre(it) = avlearnspre(it) + log(ypregrid(it,iz,ie))*zdist(it,iz)*edist(it,ie)
        avlearnspre2(it) = avlearnspre2(it) + log(ypregrid(it,iz,ie)**2)*zdist(it,iz)*edist(it,ie)
        avlearnspost(it) = avlearnspost(it) + log(ygrid(it,iz,ie))*zdist(it,iz)*edist(it,ie)
        avlearnspost2(it) = avlearnspost2(it) + log(ygrid(it,iz,ie)**2)*zdist(it,iz)*edist(it,ie)
     END DO
    END DO
END DO
varearnspre = avearnspre2 - avearnspre**2
varearnspost = avearnspost2 - avearnspost**2
varlearnspre = avlearnspre2 - avlearnspre**2
varlearnspost = avlearnspost2 - avlearnspost**2

FnTaxParam = ltottax/ltotlabincpre - targetTaxToLabinc
totlabincpre = ltotlabincpre
totlabincpost = ltotlabincpost

IF (Display==1)write(*,*) ' Tax revenue / Pre-tax labor income: ', (ltottax/ltotlabincpre)*100, '%'
    
END FUNCTION FnTaxParam

