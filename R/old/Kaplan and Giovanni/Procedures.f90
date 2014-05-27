MODULE Procedures

USE Parameters
USE Globals

IMPLICIT NONE
CONTAINS

!-----------------------------------------------

SUBROUTINE LinInterp (n,x,y,ni,xi,yi)
!this does linear interpolation of (x,y) at points xi
!requires x to be sorted in ascending order
!extrapolates out of range

INTEGER, INTENT(in)		:: n,ni
REAL(long), INTENT(in)	:: x(n),y(n),xi(ni)
REAL(long), INTENT(out)	:: yi(ni)
REAL(long), DIMENSION(ni)	:: xL,xH,yL,yH
INTEGER					:: i,j,locL(ni)


DO i = 1,ni

	LocL(i) = MAXLOC(x,1,MASK=xi(i)>x)
	
	IF (xi(i)<=x(1)) THEN
		LocL(i) = 1
	END IF

	IF (LocL(i)>=n) THEN 
		LocL(i) = n-1
	END IF


	xL(i) = x(locL(i))
	xH(i) = x(locL(i)+1)
	yL(i) = y(locL(i))
	yH(i) = y(locL(i)+1)
	
	yi(i) = yL(i) + (xi(i)-xL(i))*((yH(i)-yL(i))/(xH(i)-xL(i)))

END DO


END SUBROUTINE LinInterp

!-----------------------------------------------

SUBROUTINE LinInterp1 (n,x,y,xi,yi)
!this does linear interpolation of (x,y) at points only point,xi
!requires x to be sorted in ascending order
!extrapolates out of range

INTEGER, INTENT(in)		:: n
REAL(long), INTENT(in)	:: x(:),y(:),xi
REAL(long), INTENT(out)	:: yi
REAL(long)	            :: xL,xH,yL,yH
INTEGER					:: locL

LocL = MAXLOC(x,1,MASK=xi>x)

IF (xi<=x(1)) THEN
	LocL = 1
END IF

IF (LocL>=n) THEN 
	LocL = n-1
END IF

xL  = x(locL)
xH  = x(locL +1)
yL  = y(locL)
yH  = y(locL +1)

yi  = yL  + (xi -xL )*((yH -yL )/(xH -xL ))

END SUBROUTINE LinInterp1

!-----------------------------------------------

SUBROUTINE FindLinProb1 (xi,x,y,p)
!this takes in xi, and finds two points in either side of in in x
!and returns the indices of them y and associated probabilities p

REAL(8), INTENT(in)	    :: x(:),xi
REAL(8), INTENT(out)	:: p(2)
INTEGER, INTENT(out)    :: y(2)

INTEGER					:: locL,n

n = size(x)
LocL = MAXLOC(x,1,MASK=xi>x)

IF (xi<=x(1)) THEN
	y(1) = 1
	y(2) = 2
	p(1) = 1.0
	p(2) = 0.0
ELSEIF (LocL>=n) THEN 
	LocL = n-1
    y(1) = n-1
    y(2) = n
    p(1) = 0.0
    p(2) = 1.0
ELSE
    y(1) = LocL
    y(2) = LocL+1
    p(1) = (xi - x(LocL))/real(x(LocL+1)-x(LocL))
    p(2) = 1.0 - p(1)
END IF

END SUBROUTINE FindLinProb1

!----------------------------------------------
SUBROUTINE WriteMatrix(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
REAL(long),INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1,i2

WRITE(UNIT=lstring, FMT='(I5)') n2
lstring = '('//trim(lstring) // 'F16.6)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrix
!----------------------------------------------
SUBROUTINE WriteMatrixLong(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
REAL(long),INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1,i2

WRITE(UNIT=lstring, FMT='(I5)') n2
lstring = '('//trim(lstring) // 'F20.14)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrixLong

!----------------------------------------------
SUBROUTINE WriteMatrixInteger(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
INTEGER,INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1,i2

WRITE(UNIT=lstring, FMT='(I5)') n2
lstring = '('//trim(lstring) // 'I16)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrixInteger
!----------------------------------------------
SUBROUTINE WriteMatrixCSV(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
REAL(long),INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1,i2

WRITE(UNIT=lstring, FMT='(I5)') n2-1
lstring = '('//trim(lstring) // '(F16.6,","),F16.6)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrixCSV
!----------------------------------------------
SUBROUTINE WriteMatrixCSVInteger(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
INTEGER,INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1,i2

WRITE(UNIT=lstring, FMT='(I5)') n2-1
lstring = '('//trim(lstring) // '(I16,","),I16)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrixCSVInteger
!----------------------------------------------
SUBROUTINE RandomDiscrete(Nout,Xout,Nin,Pin)
!generates Nout random draws from the integers 1 to Nin
!using probabilities in Pin

INTEGER, INTENT(in)			:: Nout,Nin
INTEGER,INTENT(out)		:: Xout(:)
REAL(long), INTENT(in)  ::Pin(:)

INTEGER			::i1,i2
REAL(long)      :: lran(Nout)

!IF(sum(Pin) .ne. 1.0) write(*,*) 'error in RandomDiscrete: Pin doesnt sum to 1.0'

CALL RANDOM_NUMBER(lran)

Xout(:) = 0
DO i1 = 1,Nout
    IF ( lran(i1) .le. Pin(1) ) THEN
        Xout(i1) = 1
    ELSE
        i2 = 2
        DO WHILE (i2 .le. Nin)
            IF ( (lran(i1) .le. SUM(Pin(1:i2)) ).and. (lran(i1) > SUM(Pin(1:i2-1)) ) ) THEN
                Xout(i1) = i2
                i2 = Nin+1
            ELSE
                i2 = i2+1
            END IF
        END DO
    END IF
END DO

END SUBROUTINE RandomDiscrete

!--------------------------------------------------------------
SUBROUTINE RandomDiscrete1(Xout,Nin,Pin)
!generates Nout random draws from the integers 1 to Nin
!using probabilities in Pin

INTEGER, INTENT(in)			:: Nin
INTEGER,INTENT(out)		:: Xout
REAL(long), INTENT(in)  ::Pin(:)

INTEGER			::i2
REAL(long)      :: lran

!IF(sum(Pin) .ne. 1.0) write(*,*) 'error in RandomDiscrete: Pin doesnt sum to 1.0'

CALL RANDOM_NUMBER(lran)

Xout = 0
IF ( lran .le. Pin(1) ) THEN
    Xout = 1
ELSE
    i2 = 2
    DO WHILE (i2 .le. Nin)
        IF ( (lran .le. SUM(Pin(1:i2)) ).and. (lran > SUM(Pin(1:i2-1)) ) ) THEN
            Xout = i2
            i2 = Nin+1
        ELSE
            i2 = i2+1
        END IF
    END DO
END IF


END SUBROUTINE RandomDiscrete1


!--------------------------------------------------------------
RECURSIVE SUBROUTINE quick_sort(list, order)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

IMPLICIT NONE
REAL(long), DIMENSION (:), INTENT(IN OUT)  :: list
INTEGER, DIMENSION (:), INTENT(OUT)  :: order

! Local variable
INTEGER :: i

DO i = 1, SIZE(list)
  order(i) = i
END DO

CALL quick_sort_1(1, SIZE(list))

CONTAINS

RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL(long)                :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1


SUBROUTINE interchange_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL                :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort

END MODULE Procedures