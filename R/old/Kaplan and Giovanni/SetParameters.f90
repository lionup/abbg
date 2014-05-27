SUBROUTINE SetParameters

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8) :: Veta_rho1

!EARNINGS PROCESS
Veps =  0.05
Vz0  =  0.15
rho  =  1.0
Veta_rho1 =  0.01   !Veta if rho==1

IF (rho==1.0) THEN
    Veta = Veta_rho1
ELSE    
    Veta =   (1-rho**2)*(Twork*Veta_rho1 - Vz0*(rho**(2*Twork)-1) ) /(1-rho**(2*Twork))
END IF
Vetavec = Veta

!INTEREST RATE
R = 1.03      !annual gross interest rate

!PREFERENCE PARMETERS
gam = 2.0
bet =   1/R
qprefb = 120000.0 !10000000.0  !bliss point
qprefa = 1.0 !curvature in pref 

!BORROWING LIMIT: SET TO VERY LARGE NEGATIVE NO. FOR NBL
borrowlim = 0.0 !-100000000.0

!GOVERNMENT PARAMETERS
!gouveia strauss 
stax = 2.0e-4   !guess: is chosen optimally
ptax = 0.768
btax = 0.258

!other taxes and benefits
 cfloor = -100000000000.0  
pentax  = 0.0    !payroll tax   
rtax = 0.0 !tax on interest income
otax = 0.0 !tax on old age pensions, check option in Parameters
pencapfrac = 2.2 !cap on (pre-tax) earnings that contribute to pension index, as fraction of av pre-tax earnigs

Rnet = 1.0 + (1.0-rtax)*(R-1)      !annual after tax interest rate

END SUBROUTINE SetParameters