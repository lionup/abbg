SUBROUTINE rtnewtGrossInc(lnet,lxguess,lx1,lx2,xacc)
IMPLICIT NONE

REAL(8), INTENT(IN) :: lnet,lx1,lx2,xacc
REAL(8), INTENT(INOUT) :: lxguess

INTEGER, PARAMETER :: MAXIT=30
INTEGER :: j
REAL(8) :: df,dx,f,x1,x2,x

x  = lxguess
x1 = lx1
x2 = lx2


do j=1,MAXIT
	
	call FnGrossInc(x,lnet,f,df)
	dx=f/df
	x=x-dx
	lxguess = x
	if ((x1-x)*(x-x2) < 0.0) then
	    write(*,*) 'rtnewtGrossInc: values jumped out of brackets'
		return
    end if		
	if (abs(dx) < xacc) return
end do

write(*,*)('rtnewt: exceed maximum iterations')

END SUBROUTINE rtnewtGrossInc

