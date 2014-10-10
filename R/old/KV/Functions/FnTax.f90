REAL(8) FUNCTION FnTax(ly)
USE Parameters
USE Globals
IMPLICIT NONE

REAL(8), INTENT(IN) :: ly

FnTax = btax*(ly - (ly**(-ptax) + stax)**(-1.0/ptax)) + pentax*ly
    
END FUNCTION FnTax
