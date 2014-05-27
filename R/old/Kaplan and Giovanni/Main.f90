PROGRAM Main

USE Parameters
USE Globals

IF (RunOnLinux==0) filesep = '\'
IF (RunOnLinux==1) filesep = '/'

IF (MatchKYRatio==0 .and. EquilibriumR==0) THEN
    OutputDir = OutputBaseDir  
    CALL SetParameters
    CALL LoadData
    CALL Grids
    CALL Decisions
    CALL Simulate
    CALL SaveOutput
END IF

IF (MatchKYRatio==1 .and. EquilibriumR==0) THEN
    OutputDir = OutputBaseDir  
    CALL SetParameters
    CALL LoadData
    CALL Grids
    CALL rtsecBetaKY(bet,bet*1.01,0.01_8)
    CALL SaveOutput
END IF

IF (MatchKYRatio==0 .and. EquilibriumR==1) THEN
    OutputDir = OutputBaseDir  
    CALL SetParameters
    CALL LoadData
    CALL Grids
    CALL rtsecEqumR(R-1.0,(R-1.0)*1.01,0.01_8)
    CALL SaveOutput
END IF 


END PROGRAM Main