! --------------------------------------------------------------------
! MODULE  Sorting:
!    This module can sort a set of numbers.  The method used is
! usually referred to as "selection" method.
! --------------------------------------------------------------------

MODULE  MOD_SORT

USE SET_KND

IMPLICIT  NONE
PRIVATE   :: FindMinimum, Swap

CONTAINS

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

INTEGER(I4) FUNCTION  FindMinimum(x, Start, End)
IMPLICIT  NONE
REAL(R8), DIMENSION(1:), INTENT(IN) :: x
INTEGER(I4), INTENT(IN)             :: Start, End
REAL(R8)                            :: Minimum
INTEGER(I4)                         :: Location
INTEGER(I4)                         :: i

Minimum  = x(Start)          ! assume the first is the min
Location = Start             ! record its position
DO i = Start+1, End          ! start with next elements
IF (x(i) < Minimum) THEN  !   if x(i) less than the min?
Minimum  = x(i)        !      Yes, a new minimum found
Location = i                !      record its position
END IF
END DO
FindMinimum = Location            ! return the position
END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

SUBROUTINE  Swap(a, b)
IMPLICIT  NONE
REAL(R8), INTENT(INOUT) :: a, b
REAL(R8)                :: Temp

Temp = a
a    = b
b    = Temp
END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

SUBROUTINE  Sort(x, Size)
IMPLICIT  NONE
REAL(R8), DIMENSION(1:), INTENT(INOUT) :: x
INTEGER(I4), INTENT(IN)                   :: Size
INTEGER(I4)                               :: i
INTEGER(I4)                               :: Location

DO i = 1, Size-1             ! except for the last
Location = FindMinimum(x, i, Size)  ! find min from this to last
CALL  Swap(x(i), x(Location))  ! swap this and the minimum
END DO
END SUBROUTINE  Sort


END MODULE  MOD_SORT
