SUBROUTINE SaveOutput

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

write(*,*) ' Saving simulations to disk'

OutputDirSims = trim(OutputBaseDir)

OPEN(3, FILE = trim(OutputDirSims) // 'Twork.txt', STATUS = 'replace')
WRITE(3,*) Twork
CLOSE(3)
OPEN(3, FILE = trim(OutputDirSims) // 'Tret.txt', STATUS = 'replace')
WRITE(3,*) Tret
CLOSE(3)
OPEN(3, FILE = trim(OutputDirSims)// 'nsim.txt', STATUS = 'replace')
WRITE(3,*) nsim
CLOSE(3)

OPEN(3, FILE = trim(OutputDirSims)// 'rho.txt', STATUS = 'replace')
WRITE(3,*) rho
CLOSE(3)

OPEN(3, FILE = trim(OutputDirSims)// 'bet.txt', STATUS = 'replace')
WRITE(3,*) bet
CLOSE(3)

OPEN(3, FILE = trim(OutputDirSims) // 'Vetavec.txt', STATUS = 'replace')
CALL WriteMatrix(3,Twork-1,1,Vetavec)


OPEN(3, FILE = trim(OutputDirSims) // 'asim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Ttot+1,asim)

OPEN(3, FILE = trim(OutputDirSims)// 'ysim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Ttot,ysim)

OPEN(3, FILE = trim(OutputDirSims)// 'ypresim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Ttot,ypresim)

OPEN(3, FILE = trim(OutputDirSims)// 'xsim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Ttot,xsim)

OPEN(3, FILE = trim(OutputDirSims)// 'trsim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Ttot,trsim)

OPEN(3, FILE = trim(OutputDirSims) // 'csim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Ttot,csim)

OPEN(3, FILE = trim(OutputDirSims) // 'tsim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Twork,tsim)

OPEN(3, FILE = trim(OutputDirSims) // 'zsim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Twork,zsim)

OPEN(3, FILE = trim(OutputDirSims) // 'esim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Twork,esim)

OPEN(3, FILE = trim(OutputDirSims) // 'zsimI.txt', STATUS = 'replace')
CALL WriteMatrixCSVInteger(3,nsim,Twork,zsimI)

OPEN(3, FILE = trim(OutputDirSims) // 'esimI.txt', STATUS = 'replace')
CALL WriteMatrixCSVInteger(3,nsim,Twork,esimI)

OPEN(3, FILE = trim(OutputDirSims) // 'msim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Twork,msim)

OPEN(3, FILE = trim(OutputDirSims) // 'msimI.txt', STATUS = 'replace')
CALL WriteMatrixCSVInteger(3,nsim,Twork,msimI)

OPEN(3, FILE = trim(OutputDirSims) // 'yavsim.txt', STATUS = 'replace')
CALL WriteMatrixCSV(3,nsim,Twork,yavsim)

END SUBROUTINE SaveOutput