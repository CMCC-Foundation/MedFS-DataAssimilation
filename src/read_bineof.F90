SUBROUTINE READ_BINEOF(CFILEBIN,KPARS,&
           & KMT,IM,JM,KM,NREG,NEOF,EVAR,EVCR)

!... READ IN EOF IN BINARY FILE

  USE SET_KND
  USE IOUNITS, ONLY : IOUNERR, IOUNLOG, IOUNOUT,IOUNEOF

  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) ::  CFILEBIN
  INTEGER(I4),INTENT(IN)  :: KPARS,KMT,IM,JM,KM,NREG,NEOF
  REAL(R8),INTENT(OUT)  :: EVAR(NEOF,NREG), EVCR(KMT,NEOF,NREG)

  INTEGER(I4) :: K,NEC,NRG,NTOT,IPAR

  INTEGER(I4)                 :: IAUXA(10),IAUXS,IERR,IAUX0

  CHARACTER(LEN=70)            ::  CLCOMM
  CHARACTER(LEN=7)             ::  CLGRID
  CHARACTER(LEN=4)             ::  CLVERS
  CHARACTER(LEN=3)             ::  CEND
  INTEGER(I4),PARAMETER :: ICHKWD = 3141592

  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' READING BINARY EOF FILE'

  OPEN(IOUNEOF,FILE=CFILEBIN,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=IERR)
  IF(IERR.NE.0) CALL ABOR1(' CANNOT OPEN BINARY FILE '//TRIM(CFILEBIN))

! HEADERS

  READ(IOUNEOF) CLCOMM
  READ(IOUNEOF) CLGRID
  READ(IOUNEOF) CLVERS
  READ(IOUNEOF)

  WRITE(IOUNOUT,*) ' READING BIN FILE:'
  WRITE(IOUNOUT,*) ' '//TRIM(CLCOMM)

  WRITE(IOUNLOG,*) ' READING BIN FILE:'
  WRITE(IOUNLOG,*) ' '//TRIM(CLCOMM)
  WRITE(IOUNLOG,*) ' GRID: ',CLGRID,'  VERSION: ',CLVERS

  READ(IOUNEOF) IAUXA(1:7)

  WRITE(IOUNLOG,*) ' FOUND DIMENSIONS:'
  WRITE(IOUNLOG,*)  IAUXA(1:7)

  IF( IAUXA(1) .NE. KPARS .OR. &
   &  IAUXA(2) .NE. KMT   .OR. &
   &  IAUXA(3) .NE. IM    .OR. &
   &  IAUXA(4) .NE. JM    .OR. &
   &  IAUXA(5) .NE. KM    .OR. &
   &  IAUXA(6) .NE. NREG  .OR. &
   &  IAUXA(7) .NE. NEOF )  THEN

     WRITE(IOUNOUT,*) ' DIMENSIONS MISMATCH WHILE READING BINARY EOFS FILE'
     WRITE(IOUNOUT,*)
     WRITE(IOUNOUT,*) ' EXPECTED DIMENSIONS'
     WRITE(IOUNOUT,*)  KPARS,KMT,IM,JM,KM,NREG,NEOF

     CALL ABOR1('DIMENSION MISMATCH IN READ_BINEOF')

  ENDIF

  CALL FLUSH(IOUNOUT)
  CALL FLUSH(IOUNLOG)

  READ(IOUNEOF) IAUXA(1:KPARS+1)
  IF(IAUXA(KPARS+1) .NE. ICHKWD ) CALL ABOR1('PROBLEMS READING BINEOF')

  NTOT = 0
  DO IPAR=1,KPARS
     IF(IAUXA(IPAR).EQ. 3) NTOT = NTOT + KM
     IF(IAUXA(IPAR).EQ. 2) NTOT = NTOT + 1
  ENDDO

  IF( NTOT .NE. KMT ) CALL ABOR1('PROBLEMS DEFINING DIMENSIONS')

! BODY

  !... SET NO 1: EIGENVALUES
  WRITE(IOUNLOG,*) ' READING EIGENVALUES'
  READ(IOUNEOF) ((EVAR(NEC,NRG),NRG=1,NREG),NEC=1,NEOF),IAUXS
  IF(IAUXS .NE. ICHKWD ) CALL ABOR1('PROBLEMS READING BINEOF')

  !... SET NO 2: EIGENVECTORS
  WRITE(IOUNLOG,*) ' READING EIGENVECTORS'
  IAUX0 = 0
  DO K=1,KMT
     READ (IOUNEOF) &
     & ((EVCR(K,NEC,NRG),NRG=1,NREG),NEC=1,NEOF),IAUXS
     IAUX0 = IAUX0 + IAUXS
  ENDDO
  IF(IAUX0 .NE. KMT*ICHKWD ) CALL ABOR1('PROBLEMS READING BINEOF')
  READ(IOUNEOF) CEND
  IF (CEND(1:3) .NE. 'END' ) CALL ABOR1('PROBLEMS READING BINEOF')
  CLOSE(IOUNEOF)

END SUBROUTINE READ_BINEOF

!.....................................................................................

SUBROUTINE WRITE_BINEOF(CFILEBIN,CLCOMM,CLGRID,CLVERS,KPARS,KPARSL,&
           & KMT,IM,JM,KM,NREG,NEOF,EVAR,EVCR)

!... WRITE OUT IN EOF IN BINARY FILE -- FOR EXTERNAL USE

  USE SET_KND
  USE IOUNITS, ONLY : IOUNERR, IOUNLOG, IOUNOUT,IOUNEOF

  IMPLICIT NONE

  INTEGER(I4),INTENT(IN)  :: KPARS,KMT,IM,JM,KM,NREG,NEOF
  INTEGER(I4),INTENT(IN)  :: KPARSL(KPARS)
  REAL(R8),INTENT(IN)  :: EVAR(NEOF,NREG), EVCR(KMT,NEOF,NREG)

  INTEGER(I4) :: K,NEC,NRG,NTOT,IPAR

  INTEGER(I4)                 :: IAUXA(10),IAUXS,IERR

  CHARACTER(LEN=70)            ::  CLCOMM
  CHARACTER(LEN=60)            ::  CFILEBIN
  CHARACTER(LEN=7)             ::  CLGRID
  CHARACTER(LEN=4)             ::  CLVERS
  INTEGER(I4),PARAMETER :: ICHKWD = 3141592

  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' READING BINARY EOF FILE'

! HEADERS

  OPEN(IOUNEOF,FILE=CFILEBIN,FORM='UNFORMATTED')

! HEADERS
  !... WRITE COMMENTS AS SEPARATE RECORDS
  WRITE(IOUNEOF) CLCOMM
  WRITE(IOUNEOF) CLGRID
  WRITE(IOUNEOF) CLVERS
  WRITE(IOUNEOF)
  !... WRITE DIMENSIONS
  WRITE(IOUNEOF) KPARS,     KMT,    IM,    JM,    KM,    NREG,    NEOF
  ! LIST DIMENSIONALITY OF KPARS
  WRITE(IOUNEOF) KPARSL,ICHKWD
! BODY
  !... SET NO 1: EIGENVALUES
  WRITE(IOUNEOF) ((EVAR(NEC,NRG),NRG=1,NREG),NEC=1,NEOF),ICHKWD
  !... SET NO 2: EIGENVECTORS
  DO K=1,KMT
     WRITE(IOUNEOF) &
     & ((EVCR(K,NEC,NRG),NRG=1,NREG),NEC=1,NEOF),ICHKWD
  ENDDO
  WRITE(IOUNEOF) 'END'
  CLOSE(IOUNEOF)

END SUBROUTINE WRITE_BINEOF
