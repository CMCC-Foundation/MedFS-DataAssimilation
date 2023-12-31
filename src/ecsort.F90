MODULE ECSORT

!..   AUTHOR: SAMI SAARINEN, ECMWF, 10/02/98
!     FIXES : SAMI SAARINEN, ECMWF, 08/11/99 : SUB-ARRAYS GO NOW CORRECTLY (LOOK FOR ADDRDIFF)
!                                              GENUINE REAL(4) SORT "RE-HABILITATED"
!                                              SIZEOF_INT, _REAL4 & _REAL8 HARDCODED !
!             SAMI SAARINEN, ECMWF, 11/10/00 : REAL*4 VERSION INCLUDED (REAL_M)
!             SAMI SAARINEN, ECMWF, 28/11/03 : CALLS TO DR_HOOK ADDED MANUALLY (ON TOP OF CY28)
!             SAMI SAARINEN, ECMWF, 18/02/05 : 64-BIT INTEGER SORTING INTRODUCED (FOR CY30)
!             SAMI SAARINEN, ECMWF, 22/02/05 : USING GENUINE 64-BIT RSORT64() => ONE-PASS THROUGH DATA
!             SAMI SAARINEN, ECMWF, 06/07/05 : "CURRENT_METHOD" MADE OPENMP-THREAD AWARE (FOR MAX. # OF THREADS = NTHRDS)
!             SAMI SAARINEN, ECMWF, 07/07/05 : QUICK-SORT METHOD FINALLY ARRIVED (AND APPLICABLE TO MULTIKEYS, TOO)
!                                              QUICK-SORT THE DEFAULT FOR SCALAR MACHINES ("NON-VPP"), VPPS IS RADIX-SORT
!             SAMI SAARINEN, ECMWF, 03/07/07 : QUICK-SORT METHOD USES STABLE APPROACH (WASN'T GUARANTEED SO BEFORE)

!_________ LOCAL VERSION _________

!::: A.S.: MODIFIED FOR LOCAL USE : NO EXTERNAL MODULE, NO OMP SUPPORT

IMPLICIT NONE
SAVE
PRIVATE

INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: JPIB = SELECTED_INT_KIND(12)

INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300)
INTEGER, PARAMETER :: JPRM = SELECTED_REAL_KIND(6,37)


INTEGER(KIND=JPIM), PARAMETER :: NTHRDS = 32 ! ***NOTE: A HARDCODED MAX NUMBER OF THREADS !!!

INTEGER(KIND=JPIM), PARAMETER :: SIZEOF_INT   = 4
INTEGER(KIND=JPIM), PARAMETER :: SIZEOF_INT8  = 8
INTEGER(KIND=JPIM), PARAMETER :: SIZEOF_REAL4 = 4
INTEGER(KIND=JPIM), PARAMETER :: SIZEOF_REAL8 = 8

INTEGER(KIND=JPIM), PARAMETER :: MIN_METHOD = 1
INTEGER(KIND=JPIM), PARAMETER :: MAX_METHOD = 3

INTEGER(KIND=JPIM), PARAMETER :: RADIXSORT_METHOD = 1
INTEGER(KIND=JPIM), PARAMETER :: HEAPSORT_METHOD  = 2
INTEGER(KIND=JPIM), PARAMETER :: QUICKSORT_METHOD = 3

!-- SELECT SUCH METHOD FOR DEFAULT_METHOD, WHICH ALSO WORKS FOR MULTIKEY SORTS
!   VECTOR MACHINES SHOULD CHOOSE RADIXSORT_METHOD, OTHERS QUICKSORT_METHOD
#if defined(VPP) || defined(NECSX)
INTEGER(KIND=JPIM), PARAMETER :: DEFAULT_METHOD = RADIXSORT_METHOD
#else
INTEGER(KIND=JPIM), PARAMETER :: DEFAULT_METHOD = QUICKSORT_METHOD
#endif
INTEGER(KIND=JPIM)            :: CURRENT_METHOD = DEFAULT_METHOD

INTERFACE KEYSORT
MODULE PROCEDURE &
     &INT_KEYSORT_1D, INT_KEYSORT_2D, &
     &INT8_KEYSORT_1D, INT8_KEYSORT_2D, &
     &REAL8_KEYSORT_1D, REAL8_KEYSORT_2D, &
     &REAL4_KEYSORT_1D, REAL4_KEYSORT_2D
END INTERFACE

INTERFACE SORTING_METHOD
MODULE PROCEDURE INT_SORTING_METHOD, STR_SORTING_METHOD
END INTERFACE

PUBLIC :: KEYSORT
PUBLIC :: INIT_INDEX, GET_RANK
PUBLIC :: SORTING_METHOD

CONTAINS

!----------------------------
!--   PUBLIC SUBROUTINES   --
!----------------------------

SUBROUTINE INT_SORTING_METHOD(INEW, IOLD)
INTEGER(KIND=JPIM), INTENT(IN)  :: INEW
INTEGER(KIND=JPIM), INTENT(OUT) :: IOLD
INTEGER(KIND=JPIM) :: ITMP
INTEGER(KIND=JPIM) :: ITID
ITMP = INEW
IF (ITMP < MIN_METHOD .OR. ITMP > MAX_METHOD) ITMP = DEFAULT_METHOD
IOLD = CURRENT_METHOD
CURRENT_METHOD = ITMP
END SUBROUTINE INT_SORTING_METHOD


SUBROUTINE STR_SORTING_METHOD(CDNEW, IOLD)
CHARACTER(LEN=*), INTENT(IN) :: CDNEW
INTEGER(KIND=JPIM), INTENT(OUT) :: IOLD
CHARACTER(LEN=LEN(CDNEW)) CLNEW
CLNEW = CDNEW
CALL TOUPPERLOC(CLNEW)
SELECT CASE (CLNEW)
CASE ('RADIX')
  CALL SORTING_METHOD(RADIXSORT_METHOD, IOLD)
CASE ('HEAP')
  CALL SORTING_METHOD(HEAPSORT_METHOD, IOLD)
CASE ('QUICK')
  CALL SORTING_METHOD(QUICKSORT_METHOD, IOLD)
CASE ('DEFAULT')
  CALL SORTING_METHOD(DEFAULT_METHOD, IOLD)
CASE DEFAULT
  CALL SORTING_METHOD(DEFAULT_METHOD, IOLD)
END SELECT
END SUBROUTINE STR_SORTING_METHOD


SUBROUTINE INIT_INDEX(INDEX)
INTEGER(KIND=JPIM), INTENT(OUT):: INDEX(:)
INTEGER(KIND=JPIM) :: I, N
N = SIZE(INDEX)
DO I=1,N
  INDEX(I) = I
ENDDO
END SUBROUTINE INIT_INDEX


SUBROUTINE GET_RANK(INDEX, RANK)
INTEGER(KIND=JPIM), INTENT(IN) :: INDEX(:)
INTEGER(KIND=JPIM), INTENT(OUT):: RANK(:)
INTEGER(KIND=JPIM) :: I, N
N = MIN(SIZE(INDEX),SIZE(RANK))
DO I=1,N
  RANK(INDEX(I)) = I
ENDDO
END SUBROUTINE GET_RANK


SUBROUTINE INT_KEYSORT_1D(RC, A, N,METHOD, DESCENDING,INDEX, INIT)
INTEGER(KIND=JPIM), INTENT(OUT)           :: RC
INTEGER(KIND=JPIM), INTENT(INOUT)         :: A(:)
INTEGER(KIND=JPIM), INTENT(IN)            :: N
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: METHOD
LOGICAL, INTENT(IN), OPTIONAL  :: DESCENDING
INTEGER(KIND=JPIM), INTENT(INOUT), TARGET, OPTIONAL :: INDEX(:)
LOGICAL, INTENT(IN), OPTIONAL  :: INIT
! === END OF INTERFACE BLOCK ===
INTEGER(KIND=JPIM) :: AA(SIZE(A),1)
INTEGER(KIND=JPIM) :: IKEY
RC = 0
IF (N <= 0) GOTO 99
IF (SIZE(A) <= 0) GOTO 99
AA(:,1) = A(:)
IKEY = 1
IF (PRESENT(DESCENDING)) THEN
  IF (DESCENDING) IKEY = -1
ENDIF
CALL KEYSORT(RC, AA, N, KEY=IKEY, METHOD=METHOD, INDEX=INDEX, INIT=INIT)
A(:) = AA(:,1)
 99   CONTINUE
END SUBROUTINE INT_KEYSORT_1D


SUBROUTINE INT8_KEYSORT_1D(RC, A, N,METHOD, DESCENDING,INDEX, INIT)
INTEGER(KIND=JPIM), INTENT(OUT)           :: RC
INTEGER(KIND=JPIB), INTENT(INOUT)         :: A(:)
INTEGER(KIND=JPIM), INTENT(IN)            :: N
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: METHOD
LOGICAL, INTENT(IN), OPTIONAL  :: DESCENDING
INTEGER(KIND=JPIM), INTENT(INOUT), TARGET, OPTIONAL :: INDEX(:)
LOGICAL, INTENT(IN), OPTIONAL  :: INIT
! === END OF INTERFACE BLOCK ===
INTEGER(KIND=JPIB) :: AA(SIZE(A),1)
INTEGER(KIND=JPIM) :: IKEY
RC = 0
IF (N <= 0) GOTO 99
IF (SIZE(A) <= 0) GOTO 99
AA(:,1) = A(:)
IKEY = 1
IF (PRESENT(DESCENDING)) THEN
  IF (DESCENDING) IKEY = -1
ENDIF
CALL KEYSORT(RC, AA, N, KEY=IKEY, METHOD=METHOD, INDEX=INDEX, INIT=INIT)
A(:) = AA(:,1)
 99   CONTINUE
END SUBROUTINE INT8_KEYSORT_1D


SUBROUTINE REAL4_KEYSORT_1D(RC, A, N,METHOD, DESCENDING,INDEX, INIT)
INTEGER(KIND=JPIM), INTENT(OUT)           :: RC
REAL(KIND=JPRM), INTENT(INOUT)         :: A(:)
INTEGER(KIND=JPIM), INTENT(IN)            :: N
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: METHOD
LOGICAL, INTENT(IN), OPTIONAL  :: DESCENDING
INTEGER(KIND=JPIM), INTENT(INOUT), TARGET, OPTIONAL :: INDEX(:)
LOGICAL, INTENT(IN), OPTIONAL  :: INIT
! === END OF INTERFACE BLOCK ===
REAL(KIND=JPRM) :: AA(SIZE(A),1)
INTEGER(KIND=JPIM) :: IKEY
RC = 0
IF (N <= 0) GOTO 99
IF (SIZE(A) <= 0) GOTO 99
AA(:,1) = A(:)
IKEY = 1
IF (PRESENT(DESCENDING)) THEN
  IF (DESCENDING) IKEY = -1
ENDIF
CALL KEYSORT(RC, AA, N, KEY=IKEY, METHOD=METHOD, INDEX=INDEX, INIT=INIT)
A(:) = AA(:,1)
 99   CONTINUE
END SUBROUTINE REAL4_KEYSORT_1D


SUBROUTINE REAL8_KEYSORT_1D(RC, A, N,METHOD, DESCENDING,INDEX, INIT)
INTEGER(KIND=JPIM), INTENT(OUT)           :: RC
REAL(KIND=JPRB), INTENT(INOUT)         :: A(:)
INTEGER(KIND=JPIM), INTENT(IN)            :: N
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: METHOD
LOGICAL, INTENT(IN), OPTIONAL  :: DESCENDING
INTEGER(KIND=JPIM), INTENT(INOUT), TARGET, OPTIONAL :: INDEX(:)
LOGICAL, INTENT(IN), OPTIONAL  :: INIT
! === END OF INTERFACE BLOCK ===
REAL(KIND=JPRB) :: AA(SIZE(A),1)
INTEGER(KIND=JPIM) :: IKEY
RC = 0
IF (N <= 0) GOTO 99
IF (SIZE(A) <= 0) GOTO 99
AA(:,1) = A(:)
IKEY = 1
IF (PRESENT(DESCENDING)) THEN
  IF (DESCENDING) IKEY = -1
ENDIF
CALL KEYSORT(RC, AA, N, KEY=IKEY, METHOD=METHOD, INDEX=INDEX, INIT=INIT)
A(:) = AA(:,1)
 99   CONTINUE
END SUBROUTINE REAL8_KEYSORT_1D


SUBROUTINE INT_KEYSORT_2D(&
     &RC, A, N,&
     &KEY, MULTIKEY, METHOD,&
     &INDEX, INIT, TRANSPOSED)

INTEGER(KIND=JPIM), INTENT(OUT)           :: RC
INTEGER(KIND=JPIM), INTENT(INOUT)         :: A(:,:)
INTEGER(KIND=JPIM), INTENT(IN)            :: N
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: KEY, METHOD
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: MULTIKEY(:)
LOGICAL, INTENT(IN), OPTIONAL  :: TRANSPOSED
INTEGER(KIND=JPIM), INTENT(INOUT), TARGET, OPTIONAL :: INDEX(:)
LOGICAL, INTENT(IN), OPTIONAL  :: INIT
! === END OF INTERFACE BLOCK ===
INTEGER(KIND=JPIM), POINTER :: IINDEX(:)
INTEGER(KIND=JPIM) :: IKEY, ISTRIDE, IMETHOD
INTEGER(KIND=JPIM) :: LDA, IPTR, I, J, SDA, IDIFF
INTEGER(KIND=JPIM), ALLOCATABLE :: DATA(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IKEYS(:)
LOGICAL IINIT, DESCENDING, LLTRANS
INTEGER(KIND=JPIM) :: ITID

RC = 0
LDA = SIZE(A, DIM=1)
SDA = SIZE(A, DIM=2)
IF (N <= 0 .OR. LDA <= 0 .OR. SDA <= 0) GOTO 99

IMETHOD = CURRENT_METHOD
IF (PRESENT(METHOD)) THEN
  IMETHOD = MIN(MAX(MIN_METHOD,METHOD),MAX_METHOD)
ENDIF

IKEY = 1
IF (PRESENT(KEY)) IKEY = KEY

IF (PRESENT(MULTIKEY)) THEN
  ALLOCATE(IKEYS(SIZE(MULTIKEY)))
  IKEYS(:) = MULTIKEY(:)
ELSE
  ALLOCATE(IKEYS(1))
  IKEYS(1) = IKEY
ENDIF

!--   ONLY THE RADIX-SORT & NOW QUICK-SORT GIVE THE RESULT WE WANT WITH MULTIPLE KEYS
IF (SIZE(IKEYS) > 1 .AND. &
    IMETHOD /= RADIXSORT_METHOD .AND. &
    IMETHOD /= QUICKSORT_METHOD) IMETHOD = DEFAULT_METHOD

IINIT = .FALSE.
IF (PRESENT(INIT)) IINIT = INIT

IF (PRESENT(INDEX)) THEN
  IINDEX => INDEX(1:N)
ELSE
  ALLOCATE(IINDEX(N))
  IINIT = .TRUE.
ENDIF

IF (IINIT) CALL INIT_INDEX(IINDEX)

ISTRIDE = 1
LLTRANS = .FALSE.
IF (PRESENT(TRANSPOSED)) LLTRANS = TRANSPOSED
IF (LLTRANS) THEN
  ISTRIDE = LDA
ELSE IF (SDA >= 2 .AND. LDA >= 1) THEN
!-- CHECK FOR PRESENCE OF SUB-ARRAY AND ADJUST LDA AUTOMATICALLY
  CALL ADDRDIFF(A(1,1),A(1,2),IDIFF)
  LDA = IDIFF/SIZEOF_INT  ! THE TRUE LEADING DIMENSION; OVERRIDES SUB-ARRAY'S ONE
ENDIF

DO J=SIZE(IKEYS),1,-1
!--   SORT BY THE LEAST SIGNIFICANT KEY FIRST
  IKEY = ABS(IKEYS(J))
  IF (IKEY == 0) CYCLE

  IF (ISTRIDE == 1) THEN
    IPTR = LDA * (IKEY - 1) + 1
  ELSE
    IPTR = IKEY
  ENDIF

  DESCENDING = (IKEYS(J) < 0)
  IF (DESCENDING) THEN
    IF (ISTRIDE == 1) THEN
      A(1:N,IKEY) = -A(1:N,IKEY)
    ELSE
      A(IKEY, 1:N) = -A(IKEY, 1:N)
    ENDIF
  ENDIF

  SELECT CASE (IMETHOD)
  CASE (RADIXSORT_METHOD)
    CALL RSORT32_FUNC(11, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
  CASE (HEAPSORT_METHOD)
    IF (ISTRIDE == 1) THEN
      CALL INT_HEAPSORT(N, A(1:N, IKEY), IINDEX, RC)
    ELSE
      CALL INT_HEAPSORT(N, A(IKEY, 1:N), IINDEX, RC)
    ENDIF
  CASE (QUICKSORT_METHOD)
    CALL ECQSORT(11, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
  END SELECT

  IF (DESCENDING) THEN
    IF (ISTRIDE == 1) THEN
      A(1:N,IKEY) = -A(1:N,IKEY)
    ELSE
      A(IKEY, 1:N) = -A(IKEY, 1:N)
    ENDIF
  ENDIF
ENDDO

DEALLOCATE(IKEYS)

IF (.NOT.PRESENT(INDEX)) THEN
  ALLOCATE(DATA(N))

  IF (ISTRIDE == 1) THEN
    DO J=1,SDA
      DO I=1,N
        DATA(I) = A(IINDEX(I),J)
      ENDDO
      DO I=1,N
        A(I,J) = DATA(I)
      ENDDO
    ENDDO
  ELSE
    DO I=1,LDA
      DO J=1,N
        DATA(J) = A(I,IINDEX(J))
      ENDDO
      DO J=1,N
        A(I,J) = DATA(J)
      ENDDO
    ENDDO
  ENDIF

  DEALLOCATE(DATA)
  DEALLOCATE(IINDEX)
ENDIF

 99   CONTINUE
END SUBROUTINE INT_KEYSORT_2D


SUBROUTINE INT8_KEYSORT_2D(&
     &RC, A, N,&
     &KEY, MULTIKEY, METHOD,&
     &INDEX, INIT, TRANSPOSED)

INTEGER(KIND=JPIM), INTENT(OUT)           :: RC
INTEGER(KIND=JPIB), INTENT(INOUT)         :: A(:,:)
INTEGER(KIND=JPIM), INTENT(IN)            :: N
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: KEY, METHOD
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: MULTIKEY(:)
LOGICAL, INTENT(IN), OPTIONAL  :: TRANSPOSED
INTEGER(KIND=JPIM), INTENT(INOUT), TARGET, OPTIONAL :: INDEX(:)
LOGICAL, INTENT(IN), OPTIONAL  :: INIT
! === END OF INTERFACE BLOCK ===
INTEGER(KIND=JPIM), POINTER :: IINDEX(:)
INTEGER(KIND=JPIM) :: IKEY, ISTRIDE, IMETHOD
INTEGER(KIND=JPIM) :: LDA, IPTR, I, J, SDA, IDIFF
INTEGER(KIND=JPIB), ALLOCATABLE :: DATA(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IKEYS(:)
LOGICAL IINIT, DESCENDING, LLTRANS
INTEGER(KIND=JPIM) :: ITID

RC = 0
LDA = SIZE(A, DIM=1)
SDA = SIZE(A, DIM=2)
IF (N <= 0 .OR. LDA <= 0 .OR. SDA <= 0) GOTO 99

IMETHOD = CURRENT_METHOD
IF (PRESENT(METHOD)) THEN
  IMETHOD = MIN(MAX(MIN_METHOD,METHOD),MAX_METHOD)
ENDIF

IKEY = 1
IF (PRESENT(KEY)) IKEY = KEY

IF (PRESENT(MULTIKEY)) THEN
  ALLOCATE(IKEYS(SIZE(MULTIKEY)))
  IKEYS(:) = MULTIKEY(:)
ELSE
  ALLOCATE(IKEYS(1))
  IKEYS(1) = IKEY
ENDIF

!--   ONLY THE RADIX-SORT & NOW QUICK-SORT GIVE THE RESULT WE WANT WITH MULTIPLE KEYS
IF (SIZE(IKEYS) > 1 .AND. &
    IMETHOD /= RADIXSORT_METHOD .AND. &
    IMETHOD /= QUICKSORT_METHOD) IMETHOD = DEFAULT_METHOD

IINIT = .FALSE.
IF (PRESENT(INIT)) IINIT = INIT

IF (PRESENT(INDEX)) THEN
  IINDEX => INDEX(1:N)
ELSE
  ALLOCATE(IINDEX(N))
  IINIT = .TRUE.
ENDIF

IF (IINIT) CALL INIT_INDEX(IINDEX)

ISTRIDE = 1
LLTRANS = .FALSE.
IF (PRESENT(TRANSPOSED)) LLTRANS = TRANSPOSED
IF (LLTRANS) THEN
  ISTRIDE = LDA
ELSE IF (SDA >= 2 .AND. LDA >= 1) THEN
!-- CHECK FOR PRESENCE OF SUB-ARRAY AND ADJUST LDA AUTOMATICALLY
  CALL ADDRDIFF(A(1,1),A(1,2),IDIFF)
  LDA = IDIFF/SIZEOF_INT8  ! THE TRUE LEADING DIMENSION; OVERRIDES SUB-ARRAY'S ONE
ENDIF

DO J=SIZE(IKEYS),1,-1
!--   SORT BY THE LEAST SIGNIFICANT KEY FIRST
  IKEY = ABS(IKEYS(J))
  IF (IKEY == 0) CYCLE

  IF (ISTRIDE == 1) THEN
    IPTR = LDA * (IKEY - 1) + 1
  ELSE
    IPTR = IKEY
  ENDIF

  DESCENDING = (IKEYS(J) < 0)
  IF (DESCENDING) THEN
    IF (ISTRIDE == 1) THEN
      A(1:N,IKEY) = -A(1:N,IKEY)
    ELSE
      A(IKEY, 1:N) = -A(IKEY, 1:N)
    ENDIF
  ENDIF

  SELECT CASE (IMETHOD)
  CASE (RADIXSORT_METHOD)
    CALL RSORT64(14, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
!    CALL RSORT32_FUNC(14, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
  CASE (HEAPSORT_METHOD)
    IF (ISTRIDE == 1) THEN
      CALL INT8_HEAPSORT(N, A(1:N, IKEY), IINDEX, RC)
    ELSE
      CALL INT8_HEAPSORT(N, A(IKEY, 1:N), IINDEX, RC)
    ENDIF
  CASE (QUICKSORT_METHOD)
    CALL ECQSORT(14, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
  END SELECT

  IF (DESCENDING) THEN
    IF (ISTRIDE == 1) THEN
      A(1:N,IKEY) = -A(1:N,IKEY)
    ELSE
      A(IKEY, 1:N) = -A(IKEY, 1:N)
    ENDIF
  ENDIF
ENDDO

DEALLOCATE(IKEYS)

IF (.NOT.PRESENT(INDEX)) THEN
  ALLOCATE(DATA(N))

  IF (ISTRIDE == 1) THEN
    DO J=1,SDA
      DO I=1,N
        DATA(I) = A(IINDEX(I),J)
      ENDDO
      DO I=1,N
        A(I,J) = DATA(I)
      ENDDO
    ENDDO
  ELSE
    DO I=1,LDA
      DO J=1,N
        DATA(J) = A(I,IINDEX(J))
      ENDDO
      DO J=1,N
        A(I,J) = DATA(J)
      ENDDO
    ENDDO
  ENDIF

  DEALLOCATE(DATA)
  DEALLOCATE(IINDEX)
ENDIF

 99   CONTINUE
END SUBROUTINE INT8_KEYSORT_2D


SUBROUTINE REAL4_KEYSORT_2D(&
     &RC, A, N,&
     &KEY, MULTIKEY, METHOD,&
     &INDEX, INIT, TRANSPOSED)

INTEGER(KIND=JPIM), INTENT(OUT)           :: RC
REAL(KIND=JPRM), INTENT(INOUT)         :: A(:,:)
INTEGER(KIND=JPIM), INTENT(IN)            :: N
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: KEY, METHOD
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: MULTIKEY(:)
LOGICAL, INTENT(IN), OPTIONAL  :: TRANSPOSED
INTEGER(KIND=JPIM), INTENT(INOUT), TARGET, OPTIONAL :: INDEX(:)
LOGICAL, INTENT(IN), OPTIONAL  :: INIT
! === END OF INTERFACE BLOCK ===
INTEGER(KIND=JPIM), POINTER :: IINDEX(:)
INTEGER(KIND=JPIM) :: IKEY, ISTRIDE, IMETHOD
INTEGER(KIND=JPIM) :: LDA, IPTR, I, J, SDA, IDIFF
REAL(KIND=JPRM), ALLOCATABLE :: DATA(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IKEYS(:)
LOGICAL IINIT, DESCENDING, LLTRANS
INTEGER(KIND=JPIM) :: ITID

RC = 0
LDA = SIZE(A, DIM=1)
SDA = SIZE(A, DIM=2)
IF (N <= 0 .OR. LDA <= 0 .OR. SDA <= 0) GOTO 99

IMETHOD = CURRENT_METHOD
IF (PRESENT(METHOD)) THEN
  IMETHOD = MIN(MAX(MIN_METHOD,METHOD),MAX_METHOD)
ENDIF

IKEY = 1
IF (PRESENT(KEY)) IKEY = KEY

IF (PRESENT(MULTIKEY)) THEN
  ALLOCATE(IKEYS(SIZE(MULTIKEY)))
  IKEYS(:) = MULTIKEY(:)
ELSE
  ALLOCATE(IKEYS(1))
  IKEYS(1) = IKEY
ENDIF

!--   ONLY THE RADIX-SORT & NOW QUICK-SORT GIVE THE RESULT WE WANT WITH MULTIPLE KEYS
IF (SIZE(IKEYS) > 1 .AND. &
    IMETHOD /= RADIXSORT_METHOD .AND. &
    IMETHOD /= QUICKSORT_METHOD) IMETHOD = DEFAULT_METHOD

IINIT = .FALSE.
IF (PRESENT(INIT)) IINIT = INIT

IF (PRESENT(INDEX)) THEN
  IINDEX => INDEX(1:N)
ELSE
  ALLOCATE(IINDEX(N))
  IINIT = .TRUE.
ENDIF

IF (IINIT) CALL INIT_INDEX(IINDEX)

ISTRIDE = 1
LLTRANS = .FALSE.
IF (PRESENT(TRANSPOSED)) LLTRANS = TRANSPOSED
IF (LLTRANS) THEN
  ISTRIDE = LDA
ELSE IF (SDA >= 2 .AND. LDA >= 1) THEN
!-- CHECK FOR PRESENCE OF SUB-ARRAY AND ADJUST LDA AUTOMATICALLY
  CALL ADDRDIFF(A(1,1),A(1,2),IDIFF)
  LDA = IDIFF/SIZEOF_REAL4  ! THE TRUE LEADING DIMENSION; OVERRIDES SUB-ARRAY'S ONE
ENDIF

DO J=SIZE(IKEYS),1,-1
!--   SORT BY LEAST SIGNIFICANT KEY FIRST
  IKEY = ABS(IKEYS(J))
  IF (IKEY == 0) CYCLE

  IF (ISTRIDE == 1) THEN
    IPTR = LDA * (IKEY - 1) + 1
  ELSE
    IPTR = IKEY
  ENDIF

  DESCENDING = (IKEYS(J) < 0)
  IF (DESCENDING) THEN
    IF (ISTRIDE == 1) THEN
      A(1:N,IKEY) = -A(1:N,IKEY)
    ELSE
      A(IKEY, 1:N) = -A(IKEY, 1:N)
    ENDIF
  ENDIF

  SELECT CASE (IMETHOD)
  CASE (RADIXSORT_METHOD)
    CALL RSORT32_FUNC(13, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
  CASE (HEAPSORT_METHOD)
    IF (ISTRIDE == 1) THEN
      CALL REAL4_HEAPSORT(N, A(1:N, IKEY), IINDEX, RC)
    ELSE
      CALL REAL4_HEAPSORT(N, A(IKEY, 1:N), IINDEX, RC)
    ENDIF
  CASE (QUICKSORT_METHOD)
    CALL ECQSORT(13, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
  END SELECT

  IF (DESCENDING) THEN
    IF (ISTRIDE == 1) THEN
      A(1:N,IKEY) = -A(1:N,IKEY)
    ELSE
      A(IKEY, 1:N) = -A(IKEY, 1:N)
    ENDIF
  ENDIF
ENDDO

DEALLOCATE(IKEYS)

IF (.NOT.PRESENT(INDEX)) THEN
  ALLOCATE(DATA(N))

  IF (ISTRIDE == 1) THEN
    DO J=1,SDA
      DO I=1,N
        DATA(I) = A(IINDEX(I),J)
      ENDDO
      DO I=1,N
        A(I,J) = DATA(I)
      ENDDO
    ENDDO
  ELSE
    DO I=1,LDA
      DO J=1,N
        DATA(J) = A(I,IINDEX(J))
      ENDDO
      DO J=1,N
        A(I,J) = DATA(J)
      ENDDO
    ENDDO
  ENDIF

  DEALLOCATE(DATA)
  DEALLOCATE(IINDEX)
ENDIF

 99   CONTINUE
END SUBROUTINE REAL4_KEYSORT_2D


SUBROUTINE REAL8_KEYSORT_2D(&
     &RC, A, N,&
     &KEY, MULTIKEY, METHOD,&
     &INDEX, INIT, TRANSPOSED)

INTEGER(KIND=JPIM), INTENT(OUT)           :: RC
REAL(KIND=JPRB), INTENT(INOUT)         :: A(:,:)
INTEGER(KIND=JPIM), INTENT(IN)            :: N
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: KEY, METHOD
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: MULTIKEY(:)
LOGICAL, INTENT(IN), OPTIONAL  :: TRANSPOSED
INTEGER(KIND=JPIM), INTENT(INOUT), TARGET, OPTIONAL :: INDEX(:)
LOGICAL, INTENT(IN), OPTIONAL  :: INIT
! === END OF INTERFACE BLOCK ===
INTEGER(KIND=JPIM), POINTER :: IINDEX(:)
INTEGER(KIND=JPIM) :: IKEY, ISTRIDE, IMETHOD
INTEGER(KIND=JPIM) :: LDA, IPTR, I, J, SDA, IDIFF
REAL(KIND=JPRB), ALLOCATABLE :: DATA(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IKEYS(:)
LOGICAL IINIT, DESCENDING, LLTRANS
INTEGER(KIND=JPIM) :: ITID

RC = 0
LDA = SIZE(A, DIM=1)
SDA = SIZE(A, DIM=2)
IF (N <= 0 .OR. LDA <= 0 .OR. SDA <= 0) GOTO 99

IMETHOD = CURRENT_METHOD
IF (PRESENT(METHOD)) THEN
  IMETHOD = MIN(MAX(MIN_METHOD,METHOD),MAX_METHOD)
ENDIF

IKEY = 1
IF (PRESENT(KEY)) IKEY = KEY

IF (PRESENT(MULTIKEY)) THEN
  ALLOCATE(IKEYS(SIZE(MULTIKEY)))
  IKEYS(:) = MULTIKEY(:)
ELSE
  ALLOCATE(IKEYS(1))
  IKEYS(1) = IKEY
ENDIF

!--   ONLY THE RADIX-SORT & NOW QUICK-SORT GIVE THE RESULT WE WANT WITH MULTIPLE KEYS
IF (SIZE(IKEYS) > 1 .AND. &
    IMETHOD /= RADIXSORT_METHOD .AND. &
    IMETHOD /= QUICKSORT_METHOD) IMETHOD = DEFAULT_METHOD

IINIT = .FALSE.
IF (PRESENT(INIT)) IINIT = INIT

IF (PRESENT(INDEX)) THEN
  IINDEX => INDEX(1:N)
ELSE
  ALLOCATE(IINDEX(N))
  IINIT = .TRUE.
ENDIF

IF (IINIT) CALL INIT_INDEX(IINDEX)

ISTRIDE = 1
LLTRANS = .FALSE.
IF (PRESENT(TRANSPOSED)) LLTRANS = TRANSPOSED
IF (LLTRANS) THEN
  ISTRIDE = LDA
ELSE IF (SDA >= 2 .AND. LDA >= 1) THEN
!-- CHECK FOR PRESENCE OF SUB-ARRAY AND ADJUST LDA AUTOMATICALLY
  CALL ADDRDIFF(A(1,1),A(1,2),IDIFF)
  LDA = IDIFF/SIZEOF_REAL8  ! THE TRUE LEADING DIMENSION; OVERRIDES SUB-ARRAY'S ONE
ENDIF

DO J=SIZE(IKEYS),1,-1
!--   SORT BY LEAST SIGNIFICANT KEY FIRST
  IKEY = ABS(IKEYS(J))
  IF (IKEY == 0) CYCLE

  IF (ISTRIDE == 1) THEN
    IPTR = LDA * (IKEY - 1) + 1
  ELSE
    IPTR = IKEY
  ENDIF

  DESCENDING = (IKEYS(J) < 0)
  IF (DESCENDING) THEN
    IF (ISTRIDE == 1) THEN
      A(1:N,IKEY) = -A(1:N,IKEY)
    ELSE
      A(IKEY, 1:N) = -A(IKEY, 1:N)
    ENDIF
  ENDIF

  SELECT CASE (IMETHOD)
  CASE (RADIXSORT_METHOD)
    CALL RSORT64(12, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
!    CALL RSORT32_FUNC(12, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
  CASE (HEAPSORT_METHOD)
    IF (ISTRIDE == 1) THEN
      CALL REAL8_HEAPSORT(N, A(1:N, IKEY), IINDEX, RC)
    ELSE
      CALL REAL8_HEAPSORT(N, A(IKEY, 1:N), IINDEX, RC)
    ENDIF
  CASE (QUICKSORT_METHOD)
    CALL ECQSORT(12, N, ISTRIDE, IPTR, A(1,1), IINDEX(1), 1, RC)
  END SELECT

  IF (DESCENDING) THEN
    IF (ISTRIDE == 1) THEN
      A(1:N,IKEY) = -A(1:N,IKEY)
    ELSE
      A(IKEY, 1:N) = -A(IKEY, 1:N)
    ENDIF
  ENDIF
ENDDO

DEALLOCATE(IKEYS)

IF (.NOT.PRESENT(INDEX)) THEN
  ALLOCATE(DATA(N))

  IF (ISTRIDE == 1) THEN
    DO J=1,SDA
      DO I=1,N
        DATA(I) = A(IINDEX(I),J)
      ENDDO
      DO I=1,N
        A(I,J) = DATA(I)
      ENDDO
    ENDDO
  ELSE
    DO I=1,LDA
      DO J=1,N
        DATA(J) = A(I,IINDEX(J))
      ENDDO
      DO J=1,N
        A(I,J) = DATA(J)
      ENDDO
    ENDDO
  ENDIF

  DEALLOCATE(DATA)
  DEALLOCATE(IINDEX)
ENDIF

 99   CONTINUE
END SUBROUTINE REAL8_KEYSORT_2D

!-----------------------------
!--   PRIVATE SUBROUTINES   --
!-----------------------------

SUBROUTINE INT_HEAPSORT(N, A, INDEX, RC)

INTEGER(KIND=JPIM), INTENT(IN)  :: N
INTEGER(KIND=JPIM), INTENT(IN)  :: A(:)
INTEGER(KIND=JPIM), INTENT(INOUT) :: INDEX(:), RC
INTEGER(KIND=JPIM) :: I,J,RIGHT,LEFT, IDX
INTEGER(KIND=JPIM) :: TMP
RC = N
IF (N <= 1) RETURN
LEFT  = N/2+1
RIGHT = N
LOOP: DO
  IF (LEFT > 1) THEN
    LEFT = LEFT - 1
    IDX  = INDEX(LEFT)
  ELSE
    IDX = INDEX(RIGHT)
    INDEX(RIGHT) = INDEX(1)
    RIGHT = RIGHT - 1
    IF (RIGHT == 1) THEN
      INDEX(1) = IDX
      EXIT LOOP
    ENDIF
  ENDIF
  TMP = A(IDX)
  I = LEFT
  J = 2*LEFT
  DO WHILE (J <= RIGHT)
    IF (J < RIGHT) THEN
      IF (A(INDEX(J)) < A(INDEX(J+1))) J = J + 1
    ENDIF
    IF (TMP < A(INDEX(J))) THEN
      INDEX(I) = INDEX(J)
      I = J
      J = 2*J
    ELSE
      J = RIGHT + 1
    ENDIF
  ENDDO
  INDEX(I) = IDX
ENDDO LOOP
END SUBROUTINE INT_HEAPSORT


SUBROUTINE INT8_HEAPSORT(N, A, INDEX, RC)

INTEGER(KIND=JPIM), INTENT(IN)  :: N
INTEGER(KIND=JPIB), INTENT(IN)  :: A(:)
INTEGER(KIND=JPIM), INTENT(INOUT) :: INDEX(:), RC
INTEGER(KIND=JPIM) :: I,J,RIGHT,LEFT, IDX
INTEGER(KIND=JPIM) :: TMP
RC = N
IF (N <= 1) RETURN
LEFT  = N/2+1
RIGHT = N
LOOP: DO
  IF (LEFT > 1) THEN
    LEFT = LEFT - 1
    IDX  = INDEX(LEFT)
  ELSE
    IDX = INDEX(RIGHT)
    INDEX(RIGHT) = INDEX(1)
    RIGHT = RIGHT - 1
    IF (RIGHT == 1) THEN
      INDEX(1) = IDX
      EXIT LOOP
    ENDIF
  ENDIF
  TMP = A(IDX)
  I = LEFT
  J = 2*LEFT
  DO WHILE (J <= RIGHT)
    IF (J < RIGHT) THEN
      IF (A(INDEX(J)) < A(INDEX(J+1))) J = J + 1
    ENDIF
    IF (TMP < A(INDEX(J))) THEN
      INDEX(I) = INDEX(J)
      I = J
      J = 2*J
    ELSE
      J = RIGHT + 1
    ENDIF
  ENDDO
  INDEX(I) = IDX
ENDDO LOOP
END SUBROUTINE INT8_HEAPSORT


SUBROUTINE REAL4_HEAPSORT(N, A, INDEX, RC)

INTEGER(KIND=JPIM), INTENT(IN)  :: N
REAL(KIND=JPRM), INTENT(IN)  :: A(:)
INTEGER(KIND=JPIM), INTENT(INOUT) :: INDEX(:), RC
INTEGER(KIND=JPIM) :: I,J,RIGHT,LEFT, IDX
REAL(KIND=JPRM) :: TMP
RC = N
IF (N <= 1) RETURN
LEFT  = N/2+1
RIGHT = N
LOOP: DO
  IF (LEFT > 1) THEN
    LEFT = LEFT - 1
    IDX  = INDEX(LEFT)
  ELSE
    IDX = INDEX(RIGHT)
    INDEX(RIGHT) = INDEX(1)
    RIGHT = RIGHT - 1
    IF (RIGHT == 1) THEN
      INDEX(1) = IDX
      EXIT LOOP
    ENDIF
  ENDIF
  TMP = A(IDX)
  I = LEFT
  J = 2*LEFT
  DO WHILE (J <= RIGHT)
    IF (J < RIGHT) THEN
      IF (A(INDEX(J)) < A(INDEX(J+1))) J = J + 1
    ENDIF
    IF (TMP < A(INDEX(J))) THEN
      INDEX(I) = INDEX(J)
      I = J
      J = 2*J
    ELSE
      J = RIGHT + 1
    ENDIF
  ENDDO
  INDEX(I) = IDX
ENDDO LOOP
END SUBROUTINE REAL4_HEAPSORT


SUBROUTINE REAL8_HEAPSORT(N, A, INDEX, RC)

INTEGER(KIND=JPIM), INTENT(IN)  :: N
REAL(KIND=JPRB), INTENT(IN)  :: A(:)
INTEGER(KIND=JPIM), INTENT(INOUT) :: INDEX(:), RC
INTEGER(KIND=JPIM) :: I,J,RIGHT,LEFT, IDX
REAL(KIND=JPRB) :: TMP
RC = N
IF (N <= 1) RETURN
LEFT  = N/2+1
RIGHT = N
LOOP: DO
  IF (LEFT > 1) THEN
    LEFT = LEFT - 1
    IDX  = INDEX(LEFT)
  ELSE
    IDX = INDEX(RIGHT)
    INDEX(RIGHT) = INDEX(1)
    RIGHT = RIGHT - 1
    IF (RIGHT == 1) THEN
      INDEX(1) = IDX
      EXIT LOOP
    ENDIF
  ENDIF
  TMP = A(IDX)
  I = LEFT
  J = 2*LEFT
  DO WHILE (J <= RIGHT)
    IF (J < RIGHT) THEN
      IF (A(INDEX(J)) < A(INDEX(J+1))) J = J + 1
    ENDIF
    IF (TMP < A(INDEX(J))) THEN
      INDEX(I) = INDEX(J)
      I = J
      J = 2*J
    ELSE
      J = RIGHT + 1
    ENDIF
  ENDDO
  INDEX(I) = IDX
ENDDO LOOP
END SUBROUTINE REAL8_HEAPSORT

SUBROUTINE TOUPPERLOC(CDS)
CHARACTER(LEN=*), INTENT(INOUT) :: CDS
INTEGER, PARAMETER :: ICH_A = ICHAR('A')
INTEGER, PARAMETER :: ICHA  = ICHAR('A')
INTEGER, PARAMETER :: ICHZ  = ICHAR('Z')
INTEGER :: I, ICH, NEW_ICH
CHARACTER(LEN=1) CH
DO I=1,LEN(CDS)
  CH = CDS(I:I)
  ICH = ICHAR(CH)
  IF ( ICH >= ICHA .AND. ICH <= ICHZ ) THEN
    NEW_ICH = ICH + (ICH_A - ICHA)
    CH = CHAR(NEW_ICH)
    CDS(I:I) = CH
  ENDIF
ENDDO
END SUBROUTINE TOUPPERLOC


END MODULE ECSORT
