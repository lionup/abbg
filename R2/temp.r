SUBROUTINE Grids

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER         :: it,ia,iz,ie,ip,im,in,lyavorder(nsim)
REAL(8)         :: lxa,lxb,lxc,lfa,lfb,lfc,lm,lfval,ltemp1,lannfact,lwidth1,lwidth2,lyavsim(nsim),lperc5,lperc50,lperc95
REAL(8), EXTERNAL        :: GOLDEN, FnGridTrans, FnGridPerm,FnTax,FnSSParam
INTEGER   :: iseed(1)

#Transitory shocks
lxa = 1.0e-8
lxb = 1.0
CALL MNBRAK(lxa,lxb,lxc,lfa,lfb,lfc,FnGridTrans)
lfval = GOLDEN(lxa,lxb,lxc,FnGridTrans,1.0e-2_8,lm)
lfval = FnGridTrans(lm)
    
#Permanent component
varz(1) = Vz0
DO it = 2, Twork
#  varz(it) = (rho**2)*varz(it-1) + Veta
  varz(it) = (rho**2)*varz(it-1) + Vetavec(it-1)
END DO  
lxa = 1.0e-8
lxb = 1.0
CALL MNBRAK(lxa,lxb,lxc,lfa,lfb,lfc,FnGridPerm)
lfval = GOLDEN(lxa,lxb,lxc,FnGridPerm,1.0e-2_8,lm)
lfval = FnGridPerm(lm)

#Earnings
IF (PreTaxIncome==0) THEN
   CALL zbrentTaxParamNet(1.0e-8_8,1.0_8,0.000001_8,stax)
   #write(*,*) ' Found required tax parameter: ', stax
ELSEIF (PreTaxIncome==1) THEN
   CALL zbrentTaxParam(1.0e-8_8,1.0_8,0.000001_8,stax)
   #write(*,*) ' Found required tax parameter: ', stax
END IF

#Mean pre-tax earnings grid for pension
pencap = pencapfrac* sum(avearnspre)/real(Twork)

#simulate to get distribution of average pre-tax incomes
iseed(1) = 6969
CALL RANDOM_SEED(put = iseed)   
CALL RandomDiscrete(nsim, zsimI(:,1), ngpz, zdist(1,:))
CALL RandomDiscrete(nsim, esimI(:,1), ngpe, edist(1,:))
DO in = 1,nsim
    ypresim(in,1) = ypregrid(1,zsimI(in,1),esimI(in,1))
    yavsim(in,1) = min(ypresim(in,1),pencap)
    DO it = 2,Twork
        CALL RandomDiscrete1(zsimI(in,it), ngpz, ztrans(it-1,zsimI(in,it-1),:))
        CALL RandomDiscrete1(esimI(in,it), ngpe, edist(it,:))
        ypresim(in,it) = ypregrid(it,zsimI(in,it),esimI(in,it))
        yavsim(in,it) = ((it-1)*yavsim(in,it-1) + min(ypresim(in,it),pencap))/real(it)
    END DO    
END DO
        

#equally spaced between 5th perc and median, and median and 95th perc
DO it = 1,Twork
    lyavsim = yavsim(:,it)
    call quick_sort(lyavsim,lyavorder)
    lperc5 = lyavsim(nsim*0.05)
    lperc50 = lyavsim(nsim*0.5)
    lperc95 = lyavsim(nsim*0.95)
    
    mgrid(it,1) = lperc5
    mgrid(it,(1+ngpm)/2) = lperc50
    mgrid(it,ngpm) = lperc95

    lwidth1 = (mgrid(it,(1+ngpm)/2) - mgrid(it,1))/real((ngpm-1)/2)
    lwidth2 = (mgrid(it,ngpm) - mgrid(it,(1+ngpm)/2))/real((ngpm-1)/2)
    
    DO im = 2,(ngpm-1)/2
        mgrid(it,im) = mgrid(it,im-1) +lwidth1
    END DO
    DO im = (ngpm+1)/2 +1,ngpm-1
        mgrid(it,im) = mgrid(it,im-1) +lwidth2
    END DO
END DO

#Pensions
CALL rtsecSSParam(0.1_8,1.5_8,1.0e-3_8)
    
#Natural borrowing limits
#Last period
DO ip = 1, ngpp
  nblret(Tret,ip) = -pgrid(Tret,ip)/Rnet + 0.1
END DO
#Retired
DO it = Tret-1,1,-1
  DO ip = 1, ngpp
        nblret(it,ip) = annprem(it)*nblret(it+1,ip)/Rnet-pgrid(it,ip)/Rnet + 0.1
    END DO
END DO
#Last Working Period
nbl(Twork) = nblret(1,pind(1,1,1))/Rnet - ygrid(Twork,1,1)/Rnet + 0.1
#Working
DO it = Twork-1,1,-1
    nbl(it) = nbl(it+1)/Rnet - ygrid(it,1,1)/Rnet +0.1
END DO


#Assets: 0 to amax - exponentially spaced grid - pexpgrid is parametr
DO it = 1,Twork
    agrid(it,1) = max(borrowlim,nbl(it))
    DO ia = 2,ngpa
      agrid(it,ia) = agrid(it,ia-1) + exp(pexpgrid*(ia-1))
    END DO
    ltemp1 = (amax-agrid(it,1))/(agrid(it,ngpa)-agrid(it,1))
    agrid(it,:) = agrid(it,1) + ltemp1*(agrid(it,:)-agrid(it,1))
END DO
DO it = 1,Tret
    DO ip = 1,ngpp
        agridret(it,1,ip) = max(borrowlim,nblret(it,ip))
        DO ia = 2,ngpa
          agridret(it,ia,ip) = agridret(it,ia-1,ip) + exp(pexpgrid*(ia-1))
        END DO
        ltemp1 = (amax-agridret(it,1,ip))/(agridret(it,ngpa,ip)-agridret(it,1,ip))
        agridret(it,:,ip) = agridret(it,1,ip) + ltemp1*(agridret(it,:,ip)-agridret(it,1,ip))
    END DO
END DO

#Cash in Hand
DO ia = 1,ngpa
    DO iz = 1,ngpz
    DO ie = 1,ngpe
        DO it = 1,Twork-1
            tran(it,ia,iz,ie) = max(cfloor - (Rnet*agrid(it,ia) + ygrid(it,iz,ie) - agrid(it+1,1)), 0.0)
            xgrid(it,ia,iz,ie) = Rnet*agrid(it,ia) + ygrid(it,iz,ie) + tran(it,ia,iz,ie)        
        END DO
        tran(Twork,ia,iz,ie) = max(cfloor - (Rnet*agrid(Twork,ia) + ygrid(Twork,iz,ie) - agridret(1,1,pind(1,iz,ie))), 0.0)
        xgrid(Twork,ia,iz,ie) = Rnet*agrid(Twork,ia) + ygrid(Twork,iz,ie) + tran(Twork,ia,iz,ie)        
        
    END DO
    END DO

    DO ip = 1,ngpp
        DO it =1,Tret-1
            tranret(it,ia,ip) = max(cfloor - (Rnet*agridret(it,ia,ip) + pgrid(it,ip) - agridret(it+1,1,ip) ),0.0 ) 
            xgridret(it,ia,ip) = Rnet*agridret(it,ia,ip) + pgrid(it,ip) +tranret(it,ia,ip)
        END DO
        tranret(Tret,ia,ip) = max(cfloor - (Rnet*agridret(Tret,ia,ip) + pgrid(Tret,ip) ),0.0 ) 
        xgridret(Tret,ia,ip) = Rnet*agridret(Tret,ia,ip) + pgrid(Tret,ip) +tranret(Tret,ia,ip)
    END DO
END DO

IF (Display==1) write(*,*) 'Finished forming grids'
END SUBROUTINE Grids