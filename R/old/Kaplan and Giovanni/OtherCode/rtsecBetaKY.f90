SUBROUTINE rtsecBetaKY(lx1,lx2,facc)

IMPLICIT NONE

REAL(8), INTENT(IN) :: lx1,lx2,facc
REAL(8) :: rtsec,x1,x2
INTEGER, PARAMETER :: MAXIT=30
INTEGER :: j
REAL(8) :: dx,f,fl,xl,ltemp

REAL(8), EXTERNAL        :: FnBetaKY

x1 = lx1
x2 = lx2

fl=FnBetaKY(x1)
if (abs(fl)<facc) return

f=FnBetaKY(x2)
if (abs(f)<facc) return

if (abs(fl) < abs(f)) then
	rtsec=x1
	xl=x2
	ltemp = fl
	fl = f
	f = ltemp
else
	xl=x1
	rtsec=x2
end if

do j=1,MAXIT
	dx=(xl-rtsec)*f/(f-fl)
	xl=rtsec
	fl=f
	rtsec=rtsec+dx
	f=FnBetaKY(rtsec)
    if (abs(fl)<facc) return
end do
	
write(*,*)('rtsec: exceed maximum iterations')

END SUBROUTINE rtsecBetaKY