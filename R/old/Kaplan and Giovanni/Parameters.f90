MODULE Parameters
IMPLICIT NONE
SAVE

integer,parameter :: long = 8

!DIRECTORIES
character(len=*), parameter	:: OutputBaseDir = "Output/AgeSpecificVarNoBorrow/"		!to write output files
character(len=*), parameter	:: InputDir = "Input/"			!to read input files

!OPTIONS
integer,parameter  ::  RunOnLinux           = 1
integer,parameter  ::  Display              = 0
integer,parameter  ::  MatchKYRatio         = 1
integer,parameter  ::  EquilibriumR         = 0
integer,parameter  ::  AnnuityMarkets       = 1
integer,parameter  ::  PreTaxIncome         = 0
integer,parameter  ::  IntitialWealthDist   = 0
integer,parameter  ::  QuadraticPref        = 0
integer,parameter  ::  TaxPensionsGS        = 1
integer,parameter  ::  AgeSpecificVariances = 1 !0

!OPTIONS TO GET BACK TO OLD VERSION
integer,parameter  ::  MatchAggPensionBen   = 0 !instead of average pre-tax earnings
integer,parameter  ::  UseFinalZPension     = 0 !instead of average pre-tax earnings (only if MatchAggPensionBen=1)
integer,parameter  ::  ScaleBendPoints      = 0 !instead of pension amounts
integer,parameter  ::  BendPointsPostTax    = 0 !instead of basing bend points on pre-tax earnings
integer,parameter  ::  UseNetIncKY          = 0 !instead of pre-tax inc when calculating K/Y
integer,parameter  ::  UseNetIncToScaleSS   = 0 !instead of pre-tax inc when calibrating pension to lab inc ratio (only if MatchAggPensionBen=1)


! GRIDS DIMENSION - STATE VARIABLES
integer,parameter :: ngpe = 19 !7 			    !transitory component
integer,parameter :: ngpz = 39 !11 			    !permanent component
integer,parameter :: ngpa = 100 !50 		    !asset points
integer,parameter :: ngpm = 19  			    !average earnings points
integer,parameter :: ngpp = ngpm*ngpz*ngpe      !pension points

!DEMOGRAPHIC PARAMETERS PARAMETERS
integer,parameter	:: Twork    = 35		!working years
integer,parameter   :: Tret     = 35        !retirement years
integer,parameter   :: Ttot     = Twork+Tret        !total years

!TARGET MOMENTS
real(8), parameter   :: targetKY = 2.5
real(8), parameter   :: targetTaxToLabinc = 0.20 !0.17 
real(8), parameter   :: targetSSBenToLabinc = 0.095 !0.14 !only used if MatchAggPensionBen==1
real(8), parameter   :: targetSSAvReplacement = 0.45 

!PARAMETERS FOR GRID CONSTRUCTION
real(8), parameter   :: pexpgrid = 0.18       !approaches linear as goes to 0, approaches L shaped as goes to Inf
real(8), parameter   :: amax  = 300000.0

!SIMULATION PARAMETERS	
integer, parameter	:: nsim = 50000 !5000

END MODULE Parameters

