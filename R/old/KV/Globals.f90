MODULE Globals
USE Parameters

IMPLICIT NONE


!GLOBALS FOR GRIDS

real(8), dimension(Twork,ngpe)              :: egrid        !transitory component
real(8), dimension(Twork,ngpz)              :: zgrid        !permanent component
real(8), dimension(Twork,ngpa)              :: agrid        !assets
real(8), dimension(Twork,ngpz,ngpe)         :: ygrid,ypregrid !earnings
real(8), dimension(Twork,ngpm)              :: mgrid        !average pre-tax earnings grid
real(8), dimension(Twork,ngpa,ngpz,ngpe)    :: xgrid        !cash on hand

real(8), dimension(Tret,ngpp)               :: pgrid,ppregrid   !pension grid
real(8), dimension(Tret,ngpa,ngpp)          :: agridret        !assets
real(8), dimension(Tret,ngpa,ngpp)          :: xgridret        !assets


real(8), dimension(Twork)	        :: nbl		! natural borrowing limits
real(8), dimension(Tret,ngpp)	    :: nblret		! natural borrowing limits

!GLOBALS TO STORE DECISION RULES 
real(8), dimension(Twork,ngpa,ngpm,ngpz,ngpe)	:: con,ass,ass1,con1,muc,muc1,emuc1,val
real(8), dimension(Twork,ngpa,ngpz,ngpe)	    :: tran
real(8), dimension(Tret,ngpa,ngpp)	        :: conret,assret,conret1,assret1,mucret,mucret1,Emucret1,valret,tranret

!GLOBALS FOR EARNINGS DISTRIBUTIONS AND TRANSITION MAPRICES
real(8), dimension(Twork-1,ngpz,ngpz)	:: ztrans		
real(8), dimension(Twork,ngpz)			:: zdist
real(8), dimension(Twork)			    :: varz,varzapprox		
real(8), dimension(Twork,ngpe)		    :: edist	
real(8), dimension(Twork)		    :: kappa,avearnspre,avearnspost,avlearnspre,avlearnspost,avearnspre2,avearnspost2,avlearnspre2,avlearnspost2
real(8), dimension(Twork)		    :: varearnspre,varearnspost,varlearnspre,varlearnspost

!GLOBALS FOR PENSION SYSTEM
integer, dimension(ngpm,ngpz,ngpe)	:: pind

!PARAMETER GLOBALS
real(8)     :: Vz0,Veps,Veta,gam,bet,R,Rnet,borrowlim,cfloor,qprefa,qprefb
real(8)     :: simKY, ptax, stax, btax,rtax, otax,ssbend1,ssbend2,pentax,rho,pencap,pencapfrac,penratio
real(8)		:: Vetavec(Twork-1),varetaage(Twork-1)


!OTHER
character(len=1)   filesep
character(len=80)  OutputDir
character(len=80)  OutputDirSims
real(8), dimension(Tret)     :: surprob, annprem
real(8), dimension(Ttot)     :: popsize
real(8), dimension(2,75)     :: initwealthdist
real(8)     ::totlabincpre, totlabincpost

!GLOBALS TO STORE SIMULATION RESULTS
REAL(8), DIMENSION(nsim, Ttot)	    :: ysim,ypresim,csim,xsim,trsim
REAL(8), DIMENSION(nsim, Twork)	    :: esim,zsim,tsim,yavsim,msim
REAL(8), DIMENSION(nsim, Ttot+1)	:: asim
INTEGER, DIMENSION(nsim, Twork)	    :: esimI,zsimI,msimI
INTEGER, DIMENSION(nsim)	        :: psimI


END MODULE Globals
