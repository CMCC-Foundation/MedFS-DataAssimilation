      SUBROUTINE MATINV(A,N,D)
      USE SET_KND
      IMPLICIT NONE
      INTEGER(I4), INTENT(IN) :: N
      REAL(R8),INTENT(INOUT) :: A(*)
      REAL(R8) :: D, BIGA, HOLD
      INTEGER(I4) :: L(N), M(N), I, J, JK, KJ, KI, JR, IK, JQ, IJ, JI, &
      & JP, IZ, K, KK, NK

! ROUTINE          DMINV   SUBROUTINE  *****  [EYRE.RETCOF]
!       
! PURPOSE          TO INVERT A DOUBLE PRECISION MATRIX
!       
! VERSION          3.00,150385,J.R.EYRE
!       
! DESCRIPTION      ROUTINE TO INVERT A DOUBLE PRECISION MATRIX.
!                  METHOD: THE STANDARD GAUSS-JORDAN METHOD IS USED.
!                  THE DETERMINANT IS ALSO CALCULATED. A DETERMINANT
!                  OF ZERO INDICATES THAT THE MATRIX IS SINGULAR.
!       
! ARGUMENTS        A      R*8 INPUT MATRIX, DESTROYED IN COMPUTATION AND
!                             REPLACED BY RESULTANT INVERSE.
!                  N      I*4 ORDER OF MATRIX A
!                  D      R*8 RESULTANT DETERMINANT
!                  L      R*8 WORK VECTOR OF LENGTH N
!                  M      R*8 WORK VECTOR OF LENGTH N
!       
!     ..................................................................
!       
!       
!        ...............................................................
!       
!        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
!        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
!        STATEMENT WHICH FOLLOWS.
!       
!       REAL A,D,BIGA,HOLD
!       
!        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
!        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
!        ROUTINE.
!       
!        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
!        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  ABS IN STATEMENT
!        10 MUST BE CHANGED TO ABS.
!       
!        ...............................................................
!       
!        SEARCH FOR LARGEST ELEMENT
!       
      D=1.0_R8
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
   10 IF(ABS(BIGA)-ABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
!       
!        INTERCHANGE ROWS
!       
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI) =HOLD
!       
!        INTERCHANGE COLUMNS
!       
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI) =HOLD
!       
!        DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
!        CONTAINED IN BIGA)
!       
   45 IF(BIGA) 48,46,48
   46 D=0.0_R8
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
!       
!        REDUCE MATRIX
!       
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
!       
!        DIVIDE ROW BY PIVOT
!       
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
!       
!        PRODUCT OF PIVOTS
!       
      D=D*BIGA
!       
!        REPLACE PIVOT BY RECIPROCAL
!       
      A(KK)=1.0_R8/BIGA
   80 CONTINUE
!       
!        FINAL ROW AND COLUMN INTERCHANGE
!       
      K=N
  100 K=(K-1)
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI) =HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI) =HOLD
      GO TO 100
  150 RETURN
      END SUBROUTINE MATINV
