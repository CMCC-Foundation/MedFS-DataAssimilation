MODULE WEAKLY

USE SET_KND
USE GRD_STR
USE TLAD_VARS
USE EOF_STR
USE IOUNITS
USE MYNETCDF
USE RANDOM_NUMBERS
USE ADTEST

IMPLICIT NONE

INTEGER(I4) :: NXLR_WW,NYLR_WW,NZLR_WW,NEOF_WW
INTEGER(I4) :: NN_WEAKLY = 0
INTEGER(I4) :: NTS_WW
INTEGER(I4) :: CVWC4_S, CVWC4_E
INTEGER(I4) :: NWD

LOGICAL     :: LSQEOFW
INTEGER(I4) :: NREGW

LOGICAL :: LL_WEAKLY=.FALSE.
LOGICAL :: LL_TEST_WEAK = .FALSE.
LOGICAL :: LL_REDUCED_GRID 
LOGICAL :: LL_CONSTERR = .FALSE.
LOGICAL :: LL_WC_NOCORR= .FALSE.
LOGICAL :: LL_SYSTERR = .FALSE.

REAL(R8) :: EDECAY = 1.0_R8

INTEGER(I4), ALLOCATABLE :: WWKI(:,:), WWKJ(:,:)
INTEGER(I4), ALLOCATABLE :: WWKI2(:,:), WWKJ2(:,:)
INTEGER(I4), ALLOCATABLE, PRIVATE :: WWKK(:  )
REAL(R8), ALLOCATABLE :: WWPI(:,:), WWPJ(:,:)
REAL(R8), ALLOCATABLE, PRIVATE :: WWPK(:  )

REAL(R8), ALLOCATABLE :: WW(:,:,:,:)
REAL(R8), ALLOCATABLE :: WW_AD(:,:,:,:)
REAL(R8), ALLOCATABLE         :: WWT(:,:,:,:)
REAL(R8), ALLOCATABLE         :: WWT2(:,:,:)
REAL(R8), ALLOCATABLE,PRIVATE :: WWTMP1(:,:,:,:)
REAL(R8), ALLOCATABLE,PRIVATE :: WWTMP2(:,:,:,:)
REAL(R8), ALLOCATABLE,PRIVATE :: WWTMP3(:,:,:)

REAL(R8), ALLOCATABLE :: TERR(:,:,:)
REAL(R8), ALLOCATABLE :: SERR(:,:,:)
REAL(R8), ALLOCATABLE :: TERR_AD(:,:,:)
REAL(R8), ALLOCATABLE :: SERR_AD(:,:,:)
REAL(R8), ALLOCATABLE :: TERR_B(:,:,:)
REAL(R8), ALLOCATABLE :: SERR_B(:,:,:)

REAL(R8), ALLOCATABLE :: ROW(:), ROW_AD(:)

TYPE(EOF_T) :: ROSW

REAL(R8) :: WRADC = 4._R8
INTEGER(I4) :: NPHOR=3
INTEGER(I4),PARAMETER :: MAXP=1000

REAL(R8), ALLOCATABLE :: ZC(:,:,:)
REAL(R8), ALLOCATABLE :: WLON(:,:), WLAT(:,:)
INTEGER(I4), ALLOCATABLE :: NP(:,:),NPI(:,:,:),NPJ(:,:,:),WREG(:,:)

CONTAINS

SUBROUTINE WEAKLY_INIT

IMPLICIT NONE

INTEGER(I4) :: JI,JJ,JI2,JJ2,JP,II2,KP
INTEGER(I4) :: I,J,K,KK,JT
REAL(R8) :: DST,DSTX,DSTY,RADC2
REAL(R8) :: ZT1, ZT2

IF( .NOT. LL_4DVAR ) LL_WEAKLY=.FALSE.
IF( .NOT. LL_WEAKLY ) RETURN

IF( LL_WEAKLY .AND. LL_SSH) THEN
  CALL ABOR1('LL_WEAKLY .AND. LL_SSH : OPTION NOT YEAT SUPPORTED')
ENDIF

CALL GETNCDIM('WC4_GRID.nc','x',NXLR_WW)
CALL GETNCDIM('WC4_GRID.nc','y',NYLR_WW)
CALL GETNCDIM('WC4_GRID.nc','z',NZLR_WW)

ROSW%EOGNX=NXLR_WW
ROSW%EOGNY=NYLR_WW
ROSW%EOGNY=NZLR_WW
ROSW%NREG =NREGW
ROSW%KMT  =NZLR_WW*2
ROSW%NEOF =NEOF_WW
ROSW%EOF_FILE = 'MODERR_EOF.nc'

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) ' /// WEAKLY CONSTRAINT 4DVAR ///'

LL_REDUCED_GRID = .FALSE.
IF( NXLR_WW .NE. GRD%IM .OR. &
    NYLR_WW .NE. GRD%JM .OR. &
    NZLR_WW .NE. GRD%KM ) THEN
    LL_REDUCED_GRID = .TRUE.
ENDIF
WRITE(IOUNLOG,*) ' REDUCED GRID :',LL_REDUCED_GRID

ALLOCATE( WLON(NXLR_WW,NYLR_WW) )
ALLOCATE( WLAT(NXLR_WW,NYLR_WW) )

IF( LL_REDUCED_GRID ) THEN
  CALL GETNCVAR('WC4_GRID.nc','lon',NXLR_WW,NYLR_WW,WLON)
  CALL GETNCVAR('WC4_GRID.nc','lat',NXLR_WW,NYLR_WW,WLAT)
ELSE
  WLON=GRD%LON
  WLAT=GRD%LAT
ENDIF

NTS_WW = NTS_TRAJ

IF( LL_CONSTERR ) NTS_WW = 1

NN_WEAKLY = ROSW%NEOF * NXLR_WW * NYLR_WW * NTS_WW
NWD = NXLR_WW* NYLR_WW* NEOF_WW* NTS_WW

WRITE(IOUNLOG,*) ' DIMENSION  : ',NXLR_WW, NYLR_WW, NZLR_WW
WRITE(IOUNLOG,*) ' TIMESTEPS  : ',NTS_WW
WRITE(IOUNLOG,*) ' EOFS       : ',ROSW%NEOF
WRITE(IOUNLOG,*) ' CTR VECT   : ',NN_WEAKLY
WRITE(IOUNLOG,*) ' PHYS. SPACE: ',NWD
WRITE(IOUNLOG,*) 
WRITE(IOUNLOG,*) ' ALLOCATION...'

IF( LL_REDUCED_GRID ) THEN

 ALLOCATE( WWKI(GRD%IM,GRD%JM) )
 ALLOCATE( WWKJ(GRD%IM,GRD%JM) )
 ALLOCATE( WWKI2(GRD%IM,GRD%JM) )
 ALLOCATE( WWKJ2(GRD%IM,GRD%JM) )
 ALLOCATE( WWPI(GRD%IM,GRD%JM) )
 ALLOCATE( WWPJ(GRD%IM,GRD%JM) )
 ALLOCATE( WWKK(GRD%KM       ) )
 ALLOCATE( WWPK(GRD%KM       ) )

 CALL GETNCVAR('WC4_MAPPING.nc','ki',GRD%IM,GRD%JM,WWKI)
 CALL GETNCVAR('WC4_MAPPING.nc','kj',GRD%IM,GRD%JM,WWKJ)
 CALL GETNCVAR('WC4_MAPPING.nc','kk',GRD%KM       ,WWKK)
 CALL GETNCVAR('WC4_MAPPING.nc','pi',GRD%IM,GRD%JM,WWPI)
 CALL GETNCVAR('WC4_MAPPING.nc','pj',GRD%IM,GRD%JM,WWPJ)
 CALL GETNCVAR('WC4_MAPPING.nc','pk',GRD%KM       ,WWPK)

 WWKI2=WWKI+1
 WHERE(WWKI2.EQ.NXLR_WW+1) WWKI2=1

 WWKJ2=WWKJ+1
 WHERE(WWKJ2.EQ.NYLR_WW+1) WWKJ2=NYLR_WW

 WHERE(WWKK+1.GT.NZLR_WW) 
   WWKK=NZLR_WW-1
    WWPK=0._R8
 ENDWHERE

ENDIF

ALLOCATE ( ROW   ( NWD ) )
ALLOCATE ( ROW_AD( NWD ) )

ALLOCATE( WW(NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW) )
ALLOCATE( WWT(NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW) )
ALLOCATE( WWT2(NXLR_WW,NYLR_WW,2*NZLR_WW) )
ALLOCATE( WW_AD(NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW) )

ALLOCATE(TERR(GRD%IM,GRD%JM,GRD%KM) )
ALLOCATE(SERR(GRD%IM,GRD%JM,GRD%KM) )

ALLOCATE(TERR_B(GRD%IM,GRD%JM,GRD%KM) )
ALLOCATE(SERR_B(GRD%IM,GRD%JM,GRD%KM) )

ALLOCATE(TERR_AD(GRD%IM,GRD%JM,GRD%KM) )
ALLOCATE(SERR_AD(GRD%IM,GRD%JM,GRD%KM) )

WW=0._R8
WW_AD=0._R8
WWT=0._R8
WWT2=0._R8
TERR=0._R8
SERR=0._R8
TERR_AD=0._R8
SERR_AD=0._R8
TERR_B=0._R8
SERR_B=0._R8

CALL GETNCVAR('MODERR_TEM_B.nc','votemerr',&
& GRD%IM,GRD%JM,GRD%KM,TERR_B(:,:,:))
CALL GETNCVAR('MODERR_SAL_B.nc','vosalerr',&
& GRD%IM,GRD%JM,GRD%KM,SERR_B(:,:,:))

! EOFS - STRUCTURE

ALLOCATE ( ROSW%EVC( ROSW%NREG, ROSW%KMT, ROSW%NEOF), &
           ROSW%EVA( ROSW%NREG,           ROSW%NEOF) )

ALLOCATE ( WREG    ( NXLR_WW, NYLR_WW ) )

WRITE(IOUNLOG,*) ' READING EOFS :', ROSW%NEOF, LSQEOFW

CALL READ_NCEOF(TRIM(ROSW%EOF_FILE),ROSW%NREG, NXLR_WW,NYLR_WW,&
& ROSW%KMT,ROSW%NEOF,ROSW%EVA,ROSW%EVC,WREG,LSQEOFW )

ALLOCATE( NP(NXLR_WW,NYLR_WW) )
ALLOCATE( NPI(NXLR_WW,NYLR_WW,MAXP) )
ALLOCATE( NPJ(NXLR_WW,NYLR_WW,MAXP) )
ALLOCATE(  ZC(NXLR_WW,NYLR_WW,MAXP) )

CALL PREP_HORI

IF( LL_TEST_WEAK ) THEN

  ALLOCATE( WWTMP1(NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW) )
  ALLOCATE( WWTMP2(NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW) )
  WW_AD=0._R8
  ROW=0._R8
  CALL RANDOM_FLD(ROW,NWD)
  WW=0._R8
  WWT = 0._R8
  CALL CTRL2ERR_TL(0)
  CALL RANDOM_FLD(WW_AD,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  WRITE(IOUNLOG,*) 'XX ',MINVAL(WW),MAXVAL(WW)
  WWT = 0._R8
  ROW_AD=0._R8
  WWTMP1 = WW_AD
  CALL CTRL2ERR_AD
  WRITE(IOUNLOG,*) 'XX ',MINVAL(WW),MAXVAL(WW)
  ZT1 = DOT_PRODUCT(ROW,ROW_AD)
  ZT2 = MDOT_PRODUCT(WW,WWTMP1,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)

  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' ADJOINT TEST FOR CTRL2ERR'
  WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
  WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
  WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
  WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
  WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
  WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*)

  CALL RANDOM_FLD(ROW,NWD)
  CALL RANDOM_FLD(WW_AD,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  WW=0._R8
  ROW_AD=0._R8
  WRITE(IOUNLOG,*) 'WB',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
  WRITE(IOUNLOG,*) 'WB',MINVAL(ROW_AD),MAXVAL(ROW_AD),MINVAL(WW_AD),MAXVAL(WW_AD) 
  CALL VEOF_WW
  CALL VEOF_WW_AD
  WRITE(IOUNLOG,*) 'WA',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
  WRITE(IOUNLOG,*) 'WA',MINVAL(ROW_AD),MAXVAL(ROW_AD),MINVAL(WW_AD),MAXVAL(WW_AD) 
  ZT1 = DOT_PRODUCT(ROW,ROW_AD)
  ZT2 = MDOT_PRODUCT(WW,WW_AD,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' ADJOINT TEST FOR VEOF_WW'
  WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
  WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
  WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
  WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
  WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
  WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*)


  CALL RANDOM_FLD(WW,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  CALL RANDOM_FLD(WW_AD,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  WWTMP1= WW
  WWTMP2= WW_AD
  WRITE(IOUNLOG,*) 'WB',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
  WRITE(IOUNLOG,*) 'WB',MINVAL(ROW_AD),MAXVAL(ROW_AD),MINVAL(WW_AD),MAXVAL(WW_AD) 
  CALL HORI_WW
  CALL HORI_WW_AD
  WRITE(IOUNLOG,*) 'WA',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
  WRITE(IOUNLOG,*) 'WA',MINVAL(ROW_AD),MAXVAL(ROW_AD),MINVAL(WW_AD),MAXVAL(WW_AD) 
  ZT1 = MDOT_PRODUCT(WWTMP1,WW_AD,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  ZT2 = MDOT_PRODUCT(WW,WWTMP2,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' ADJOINT TEST FOR HORI_WW'
  WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
  WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
  WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
  WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
  WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
  WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*)

  IF(.NOT.LL_CONSTERR ) THEN
  CALL RANDOM_FLD(WW,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  CALL RANDOM_FLD(WW_AD,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  WWTMP1= WW
  WWTMP2= WW_AD
  WRITE(IOUNLOG,*) 'WB',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
  WRITE(IOUNLOG,*) 'WB',MINVAL(ROW_AD),MAXVAL(ROW_AD),MINVAL(WW_AD),MAXVAL(WW_AD) 
  CALL TIME_WW
  CALL TIME_WW_AD
  WRITE(IOUNLOG,*) 'WA',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
  WRITE(IOUNLOG,*) 'WA',MINVAL(ROW_AD),MAXVAL(ROW_AD),MINVAL(WW_AD),MAXVAL(WW_AD) 
  ZT1 = MDOT_PRODUCT(WWTMP1,WW_AD,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  ZT2 = MDOT_PRODUCT(WW,WWTMP2,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' ADJOINT TEST FOR TIME_WW'
  WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
  WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
  WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
  WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
  WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
  WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*)
  ENDIF

  DEALLOCATE( WWTMP1, WWTMP2 )

  JT=5
  TERR=0._R8
  SERR=0._R8
  WW_AD=0._R8
  WW=0._R8
  
  ALLOCATE( WWTMP1(GRD%IM*GRD%JM*2*GRD%KM,1,1,1) )
  ALLOCATE( WWTMP2(GRD%IM*GRD%JM*2*GRD%KM,1,1,1) )
  
  CALL RANDOM_FLD(WW,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  CALL RANDOM_FLD(WWTMP1(:,1,1,1),GRD%IM*GRD%JM*2*GRD%KM) 
  KK=0
  DO K=1,GRD%KM
   DO J=1,GRD%JM
    DO I=1,GRD%IM
      KK=KK+1
      TERR_AD(I,J,K)= WWTMP1(KK,1,1,1)
    ENDDO
   ENDDO
  ENDDO
  DO K=1,GRD%KM
   DO J=1,GRD%JM
    DO I=1,GRD%IM
      KK=KK+1
      SERR_AD(I,J,K)= WWTMP1(KK,1,1,1)
    ENDDO
   ENDDO
  ENDDO
  CALL  LOW2HI_TL(JT)
  CALL  LOW2HI_AD(JT)
  KK=0
  DO K=1,GRD%KM
   DO J=1,GRD%JM
    DO I=1,GRD%IM
      KK=KK+1
      WWTMP2(KK,1,1,1) = TERR(I,J,K)
    ENDDO
   ENDDO
  ENDDO
  DO K=1,GRD%KM
   DO J=1,GRD%JM
    DO I=1,GRD%IM
      KK=KK+1
      WWTMP2(KK,1,1,1) = SERR(I,J,K)
    ENDDO
   ENDDO
  ENDDO
  ZT1 = MDOT_PRODUCT(WW,WW_AD,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW)
  ZT2 =  DOT_PRODUCT(WWTMP2(:,1,1,1),WWTMP1(:,1,1,1) )

  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' ADJOINT TEST FOR LOW2HI '
  WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
  WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
  WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
  WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
  WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
  WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*)

  IF( LL_REDUCED_GRID ) THEN
  TERR=0._R8
  SERR=0._R8
  WWT2=0._R8
  ALLOCATE( WWTMP3(NXLR_WW,NYLR_WW,2*NZLR_WW) )
  CALL RANDOM_FLD(WWT2,NXLR_WW,NYLR_WW,2*NZLR_WW)
  CALL RANDOM_FLD(WWTMP1(:,1,1,1),GRD%IM*GRD%JM*2*GRD%KM) 
  WRITE(IOUNLOG,*) 'MINMAX1 : ',MINVAL( WWT2 ), MAXVAL( WWT2 )
  CALL SMPL_TL
  WWTMP3 = WWT2
  KK=0
  DO K=1,GRD%KM
   DO J=1,GRD%JM
    DO I=1,GRD%IM
      KK=KK+1
      TERR_AD(I,J,K)= WWTMP1(KK,1,1,1)
      WWTMP2(KK,1,1,1) = TERR(I,J,K)
      IF( ABS(TERR(I,J,K)) .GT. 10._R8 ) THEN
          WRITE(9191,*) 'T',I,J,K,TERR(I,J,K),WWKI(I,J),WWKJ(I,J),WWKK(K),&
          & WWPI(I,J),WWPJ(I,J),WWPK(K)
      ENDIF
    ENDDO
   ENDDO
  ENDDO
  DO K=1,GRD%KM
   DO J=1,GRD%JM
    DO I=1,GRD%IM
      KK=KK+1
      SERR_AD(I,J,K)= WWTMP1(KK,1,1,1)
      WWTMP2(KK,1,1,1) = SERR(I,J,K)
      IF( ABS(SERR(I,J,K)) .GT. 10._R8 ) THEN
          WRITE(9191,*) 'S',I,J,K,SERR(I,J,K),WWKI(I,J),WWKJ(I,J),WWKK(K),&
          & WWPI(I,J),WWPJ(I,J),WWPK(K)
      ENDIF
    ENDDO
   ENDDO
  ENDDO
  WRITE(IOUNLOG,*) 'MINMAX2 : ',MINVAL( WWTMP2 ), MAXVAL( WWTMP2 )
  WWT2 = 0._R8
  WRITE(IOUNLOG,*) 'MINMAX3 : ',MINVAL( WWTMP1 ), MAXVAL( WWTMP1 )
  CALL SMPL_AD
  WRITE(IOUNLOG,*) 'MINMAX4 : ',MINVAL( WWT2 ), MAXVAL( WWT2 )
  ZT1 = MDOT_PRODUCT(WWTMP3,WWT2,NXLR_WW,NYLR_WW,2*NZLR_WW)
  ZT2 =  DOT_PRODUCT(WWTMP2(:,1,1,1),WWTMP1(:,1,1,1) )

  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' ADJOINT TEST FOR SMPL '
  WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
  WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
  WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
  WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
  WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
  WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*)

  ENDIF

ENDIF

WRITE(IOUNLOG,*) ' WEAKLY_INIT ENDING'

END SUBROUTINE WEAKLY_INIT

SUBROUTINE PREP_HORI

IMPLICIT NONE

REAL(R8) :: RADC2,DSTX,DSTY,DST 
INTEGER(I4) :: JI,JJ,KP,JI2,JJ2,II2
REAL(R8), PARAMETER :: ZPI=3.1415926535_R8

RADC2=WRADC**2
NP=0
ZC=0._R8

WRITE(IOUNLOG,*) ' HORIZ. OPER. :', RADC2, NPHOR

DO JJ=1,NYLR_WW
   DO JI=1,NXLR_WW
     KP=0
     DO JJ2=JJ-NPHOR,JJ+NPHOR
        DO JI2=JI-NPHOR,JI+NPHOR
           IF( JJ2.LE. NYLR_WW .AND. JJ2.GE.1 ) THEN
             II2=JI2
             IF(II2.GT.NXLR_WW) II2=II2-NXLR_WW
             IF(II2.LT.1      ) II2=II2+NXLR_WW
             KP=KP+1
             NPI(JI,JJ,KP)=II2
             NPJ(JI,JJ,KP)=JJ2
             DSTX= ABS(WLON(II2,JJ2)-WLON(JI,JJ))
             IF(DSTX.GT.180._R8) DSTX=DSTX-360._R8
             DSTY= WLAT(II2,JJ2)-WLAT(JI,JJ)
             DST = ( DSTX*COS(ZPI*180._R8*(WLAT(II2,JJ2)+WLAT(JI,JJ))/2) )**2 + DSTY**2
             ZC(JI,JJ,KP) = EXP(-MIN(32._R8,DST/RADC2))
           ENDIF
        ENDDO
     ENDDO
     NP(JI,JJ)=KP
     IF(KP.EQ.0) CALL ABOR1('WEAKLY_INIT: ZERO KP')
     ZC(JI,JJ,1:KP) = ZC(JI,JJ,1:KP) / SUM( ZC(JI,JJ,1:KP) )
   ENDDO
ENDDO

END SUBROUTINE PREP_HORI

SUBROUTINE CTRL2ERR_TL(IMODE)

IMPLICIT NONE

  INTEGER(I4), INTENT(IN) :: IMODE

  WRITE(IOUNLOG,*) 'W1',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
  CALL FLUSH(IOUNLOG)
  CALL VEOF_WW
  WRITE(IOUNLOG,*) 'W2',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
  CALL FLUSH(IOUNLOG)
  IF(.NOT.LL_WC_NOCORR) THEN
    CALL HORI_WW
    WRITE(IOUNLOG,*) 'W3',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
    CALL FLUSH(IOUNLOG)
    IF(.NOT.LL_CONSTERR) CALL TIME_WW
    WRITE(IOUNLOG,*) 'W4',MINVAL(ROW),MAXVAL(ROW),MINVAL(WW),MAXVAL(WW) 
    CALL FLUSH(IOUNLOG)
  ENDIF

  IF(IMODE.EQ.1) THEN
    CALL WRITE_VARNCDF('MODERR_TEM.nc',NXLR_WW,NYLR_WW,NZLR_WW,NTS_WW,&
    & WW(:,:,1:NZLR_WW,:),'votemerr','degC','Temperature model error')
    CALL WRITE_VARNCDF('MODERR_SAL.nc',NXLR_WW,NYLR_WW,NZLR_WW,NTS_WW,&
    & WW(:,:,(NZLR_WW+1):(2*NZLR_WW),:),'vosalerr','psu','Salinity model error')
  ENDIF

END SUBROUTINE CTRL2ERR_TL

SUBROUTINE CTRL2ERR_AD

IMPLICIT NONE

  IF(.NOT.LL_WC_NOCORR) THEN
    WRITE(IOUNLOG,*) 'W1A',MINVAL(WW_AD),MAXVAL(WW_AD),MINVAL(ROW_AD),MAXVAL(ROW_AD); CALL FLUSH(IOUNLOG)
    IF(.NOT.LL_CONSTERR) CALL TIME_WW_AD
    WRITE(IOUNLOG,*) 'W2A',MINVAL(WW_AD),MAXVAL(WW_AD),MINVAL(ROW_AD),MAXVAL(ROW_AD); CALL FLUSH(IOUNLOG)
    CALL HORI_WW_AD
  ENDIF
  WRITE(IOUNLOG,*) 'W3A',MINVAL(WW_AD),MAXVAL(WW_AD),MINVAL(ROW_AD),MAXVAL(ROW_AD); CALL FLUSH(IOUNLOG)
  CALL VEOF_WW_AD
  WRITE(IOUNLOG,*) 'W4A',MINVAL(WW_AD),MAXVAL(WW_AD),MINVAL(ROW_AD),MAXVAL(ROW_AD); CALL FLUSH(IOUNLOG)

END SUBROUTINE CTRL2ERR_AD

SUBROUTINE VEOF_WW

IMPLICIT NONE

INTEGER(I4) :: JI,JJ,JK,JT,K,JE

WW = 0._R8
K=0
DO JT=1,NTS_WW
 DO JJ=1,NYLR_WW
  DO JI=1,NXLR_WW
   DO JE=1,NEOF_WW
    K=K+1
    DO JK=1,2*NZLR_WW
     WW(JI,JJ,JK,JT) = WW(JI,JJ,JK,JT) + &
     & ROW(K)*ROSW%EVC(WREG(JI,JJ),JK,JE)*ROSW%EVA(WREG(JI,JJ),JE)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE VEOF_WW

SUBROUTINE VEOF_WW_AD

IMPLICIT NONE

INTEGER(I4) :: JI,JJ,JK,JT,K,JE

ROW_AD = 0._R8
K=0
DO JT=1,NTS_WW
 DO JJ=1,NYLR_WW
  DO JI=1,NXLR_WW
   DO JE=1,NEOF_WW
    K=K+1
    DO JK=1,2*NZLR_WW
     ROW_AD(K)= ROW_AD(K) + &
     & WW_AD(JI,JJ,JK,JT)*ROSW%EVC(WREG(JI,JJ),JK,JE)*ROSW%EVA(WREG(JI,JJ),JE)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE VEOF_WW_AD

SUBROUTINE HORI_WW

IMPLICIT NONE

INTEGER(I4) :: JI,JJ,JK,JT,JP

WWT=WW
WW=0._R8
DO JT=1,NTS_WW
 DO JK=1,2*NZLR_WW
   DO JJ=1,NYLR_WW
    DO JI=1,NXLR_WW
     DO JP=1,NP(JI,JJ)
       WW(JI,JJ,JK,JT) = WW(JI,JJ,JK,JT)+ WWT(NPI(JI,JJ,JP),NPJ(JI,JJ,JP),JK,JT)*ZC(JI,JJ,JP)
     ENDDO
    ENDDO
   ENDDO
 ENDDO
ENDDO

END SUBROUTINE HORI_WW

SUBROUTINE HORI_WW_AD

IMPLICIT NONE

INTEGER(I4) :: JI,JJ,JK,JT,JP

WWT = WW_AD
WW_AD=0._R8
DO JT=NTS_WW,1,-1
 DO JK=2*NZLR_WW,1,-1
  DO JJ=NYLR_WW,1,-1
   DO JI=NXLR_WW,1,-1
    DO JP=NP(JI,JJ),1,-1
       WW_AD(NPI(JI,JJ,JP),NPJ(JI,JJ,JP),JK,JT) = WW_AD(NPI(JI,JJ,JP),NPJ(JI,JJ,JP),JK,JT)+ WWT(JI,JJ,JK,JT)*ZC(JI,JJ,JP)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE HORI_WW_AD

SUBROUTINE TIME_WW
IMPLICIT NONE

INTEGER(I4) :: JT
REAL(R8) :: ZC0,ZC1,ZCTA,ZCTB

ZC0=1._R8
ZC1=EXP(-EDECAY)
ZCTA=2._R8*ZC1+ZC0
ZCTB=1._R8*ZC1+ZC0

WWT = WW
WW = 0._R8

DO JT=1,NTS_WW
  IF(JT.EQ.1) THEN
   WW(:,:,:,JT) = ( WWT(:,:,:,JT)*ZC0+WWT(:,:,:,JT+1)*ZC1 ) / ZCTB
  ELSEIF(JT.EQ.NTS_WW) THEN
   WW(:,:,:,JT) = ( WWT(:,:,:,JT-1)*ZC1+WWT(:,:,:,JT)*ZC0) / ZCTB
  ELSE
   WW(:,:,:,JT) = ( WWT(:,:,:,JT-1)*ZC1+WWT(:,:,:,JT)*ZC0+WWT(:,:,:,JT+1)*ZC1 ) / ZCTA
  ENDIF
ENDDO

END SUBROUTINE TIME_WW

SUBROUTINE TIME_WW2
IMPLICIT NONE

INTEGER(I4) :: JT
REAL(R8) :: ZC0,ZC1,ZC2,ZCTA,ZCTB,ZCTC

ZC0=1._R8
ZC1=EXP(-EDECAY)
ZC2=EXP(-2*EDECAY)
ZCTA=2._R8*ZC1+ZC0+2*ZC2
ZCTB=1._R8*ZC1+ZC0+ZC2
ZCTC=2._R8*ZC1+ZC0+ZC2

WWT = WW
WW = 0._R8

DO JT=1,NTS_WW
  IF(JT.EQ.1) THEN
   WW(:,:,:,JT) = ( WWT(:,:,:,JT)*ZC0+WWT(:,:,:,JT+1)*ZC1+WWT(:,:,:,JT+2)*ZC2 ) / ZCTB
  ELSEIF(JT.EQ.2) THEN
   WW(:,:,:,JT) = ( WWT(:,:,:,JT-1)*ZC1+WWT(:,:,:,JT)*ZC0+WWT(:,:,:,JT+1)*ZC1+WWT(:,:,:,JT+2)*ZC2 ) / ZCTC
  ELSEIF(JT.EQ.NTS_WW) THEN
   WW(:,:,:,JT) = ( WWT(:,:,:,JT-2)*ZC2+WWT(:,:,:,JT-1)*ZC1+WWT(:,:,:,JT)*ZC0) / ZCTB
  ELSEIF(JT.EQ.NTS_WW-1) THEN
   WW(:,:,:,JT) = ( WWT(:,:,:,JT-2)*ZC2+WWT(:,:,:,JT-1)*ZC1+WWT(:,:,:,JT)*ZC0+WWT(:,:,:,JT+1)*ZC1) / ZCTB
  ELSE
   WW(:,:,:,JT) = ( WWT(:,:,:,JT-2)*ZC2+WWT(:,:,:,JT-1)*ZC1+WWT(:,:,:,JT)*ZC0+&
                & WWT(:,:,:,JT+1)*ZC1+WWT(:,:,:,JT+2)*ZC2 ) / ZCTA
  ENDIF
ENDDO

END SUBROUTINE TIME_WW2

SUBROUTINE TIME_WW_AD2
IMPLICIT NONE

INTEGER(I4) :: JT
REAL(R8) :: ZC0,ZC1,ZC2,ZCTA,ZCTB,ZCTC

ZC0=1._R8
ZC1=EXP(-EDECAY)
ZC2=EXP(-2*EDECAY)
ZCTA=2._R8*ZC1+ZC0+2*ZC2
ZCTB=1._R8*ZC1+ZC0+ZC2
ZCTC=2._R8*ZC1+ZC0+ZC2

WWT = WW_AD
WW_AD = 0._R8
DO JT=NTS_WW,1,-1
  IF(JT.EQ.1) THEN
   WW_AD(:,:,:,JT+2) = WW_AD(:,:,:,JT+2) + WWT(:,:,:,JT)*ZC2/ZCTB
   WW_AD(:,:,:,JT+1) = WW_AD(:,:,:,JT+1) + WWT(:,:,:,JT)*ZC1/ZCTB
   WW_AD(:,:,:,JT  ) = WW_AD(:,:,:,JT  ) + WWT(:,:,:,JT)*ZC0/ZCTB
  ELSEIF(JT.EQ.2) THEN
   WW_AD(:,:,:,JT+2) = WW_AD(:,:,:,JT+2) + WWT(:,:,:,JT)*ZC2/ZCTC
   WW_AD(:,:,:,JT+1) = WW_AD(:,:,:,JT+1) + WWT(:,:,:,JT)*ZC1/ZCTC
   WW_AD(:,:,:,JT  ) = WW_AD(:,:,:,JT  ) + WWT(:,:,:,JT)*ZC0/ZCTC
   WW_AD(:,:,:,JT-1) = WW_AD(:,:,:,JT-1) + WWT(:,:,:,JT)*ZC1/ZCTC
  ELSEIF(JT.EQ.NTS_WW) THEN
   WW_AD(:,:,:,JT-2) = WW_AD(:,:,:,JT-2) + WWT(:,:,:,JT)*ZC2/ZCTB
   WW_AD(:,:,:,JT-1) = WW_AD(:,:,:,JT-1) + WWT(:,:,:,JT)*ZC1/ZCTB
   WW_AD(:,:,:,JT  ) = WW_AD(:,:,:,JT  ) + WWT(:,:,:,JT)*ZC0/ZCTB
  ELSEIF(JT.EQ.NTS_WW-1) THEN
   WW_AD(:,:,:,JT-2) = WW_AD(:,:,:,JT-2) + WWT(:,:,:,JT)*ZC2/ZCTC
   WW_AD(:,:,:,JT-1) = WW_AD(:,:,:,JT-1) + WWT(:,:,:,JT)*ZC1/ZCTC
   WW_AD(:,:,:,JT  ) = WW_AD(:,:,:,JT  ) + WWT(:,:,:,JT)*ZC0/ZCTC
   WW_AD(:,:,:,JT+1) = WW_AD(:,:,:,JT+1) + WWT(:,:,:,JT)*ZC1/ZCTC
  ELSE
   WW_AD(:,:,:,JT-2) = WW_AD(:,:,:,JT-2) + WWT(:,:,:,JT)*ZC2/ZCTA
   WW_AD(:,:,:,JT-1) = WW_AD(:,:,:,JT-1) + WWT(:,:,:,JT)*ZC1/ZCTA
   WW_AD(:,:,:,JT  ) = WW_AD(:,:,:,JT  ) + WWT(:,:,:,JT)*ZC0/ZCTA
   WW_AD(:,:,:,JT+1) = WW_AD(:,:,:,JT+1) + WWT(:,:,:,JT)*ZC1/ZCTA
   WW_AD(:,:,:,JT+2) = WW_AD(:,:,:,JT+2) + WWT(:,:,:,JT)*ZC2/ZCTA
  ENDIF
ENDDO

END SUBROUTINE TIME_WW_AD2


SUBROUTINE TIME_WW_AD
IMPLICIT NONE

INTEGER(I4) :: JT
REAL(R8) :: ZC0,ZC1,ZCTA,ZCTB

ZC0=1._R8
ZC1=EXP(-EDECAY)
ZCTA=2._R8*ZC1+ZC0
ZCTB=1._R8*ZC1+ZC0

WWT = WW_AD
WW_AD = 0._R8
DO JT=NTS_WW,1,-1
  IF(JT.EQ.1) THEN
   WW_AD(:,:,:,JT+1) = WW_AD(:,:,:,JT+1) + WWT(:,:,:,JT)*ZC1/ZCTB
   WW_AD(:,:,:,JT  ) = WW_AD(:,:,:,JT  ) + WWT(:,:,:,JT)*ZC0/ZCTB
  ELSEIF(JT.EQ.NTS_WW) THEN
   WW_AD(:,:,:,JT-1) = WW_AD(:,:,:,JT-1) + WWT(:,:,:,JT)*ZC1/ZCTB
   WW_AD(:,:,:,JT  ) = WW_AD(:,:,:,JT  ) + WWT(:,:,:,JT)*ZC0/ZCTB
  ELSE
   WW_AD(:,:,:,JT-1) = WW_AD(:,:,:,JT-1) + WWT(:,:,:,JT)*ZC1/ZCTA
   WW_AD(:,:,:,JT  ) = WW_AD(:,:,:,JT  ) + WWT(:,:,:,JT)*ZC0/ZCTA
   WW_AD(:,:,:,JT+1) = WW_AD(:,:,:,JT+1) + WWT(:,:,:,JT)*ZC1/ZCTA
  ENDIF
ENDDO

END SUBROUTINE TIME_WW_AD

SUBROUTINE LOW2HI_TL(JT)
IMPLICIT NONE
INTEGER(I4), INTENT(IN) :: JT
INTEGER(I4) :: KT1, KT2
REAL(R8)    :: RT1
KT1 = KBG1(JT)
KT2 = KBG2(JT)
RT1 = RBG1(JT)

IF( LL_CONSTERR ) THEN
   WWT2 = WW(:,:,:,1)
 ELSE
   WWT2 = RT1*WW(:,:,:,KT1) + (1._R8-RT1)*WW(:,:,:,KT2)
ENDIF

IF( LL_REDUCED_GRID ) THEN
  CALL SMPL_TL
ELSE
  TERR = WWT2(:,:,1:GRD%KM)
  SERR = WWT2(:,:,(GRD%KM+1):(2*GRD%KM))
ENDIF

IF( LL_SYSTERR ) THEN
  TERR = TERR + TERR_B
  SERR = SERR + SERR_B
ENDIF

END SUBROUTINE LOW2HI_TL

SUBROUTINE LOW2HI_AD(JT)
IMPLICIT NONE
INTEGER(I4), INTENT(IN) :: JT
INTEGER(I4) :: KT1, KT2
REAL(R8)    :: RT1

KT1 = KBG1(JT)
KT2 = KBG2(JT)
RT1 = RBG1(JT)

WWT2 = 0._R8

IF( LL_REDUCED_GRID ) THEN
  CALL SMPL_AD
ELSE
  WWT2(:,:,1:GRD%KM) = TERR_AD
  WWT2(:,:,(GRD%KM+1):(2*GRD%KM)) = SERR_AD
ENDIF

WRITE(IOUNLOG,*) 'WTA',MINVAL(WWT2),MAXVAL(WWT2); CALL FLUSH(IOUNLOG)

IF( LL_CONSTERR ) THEN
  WW_AD(:,:,:,1) = WW_AD(:,:,:,1) + WWT2
ELSE
  WW_AD(:,:,:,KT1) = WW_AD(:,:,:,KT1) + RT1*WWT2
  WW_AD(:,:,:,KT2) = WW_AD(:,:,:,KT2) + (1._R8-RT1)*WWT2
ENDIF

WRITE(IOUNLOG,*) 'WTA',MINVAL(WW_AD),MAXVAL(WW_AD); CALL FLUSH(IOUNLOG)

END SUBROUTINE LOW2HI_AD

SUBROUTINE SMPL_TL

IMPLICIT NONE

INTEGER(I4) :: I,J,K

DO K=1,GRD%KM
  DO J=1,GRD%JM
    DO I=1,GRD%IM

 TERR(I,J,K) =        WWPI(I,J) *     WWPJ(I,J) *    WWPK(K) *WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)  ) + &
                  (1.-WWPI(I,J))*     WWPJ(I,J) *    WWPK(K) *WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)  ) + &
                      WWPI(I,J) * (1.-WWPJ(I,J))*    WWPK(K) *WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)  ) + &
                  (1.-WWPI(I,J))* (1.-WWPJ(I,J))*    WWPK(K) *WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)  ) + &
                      WWPI(I,J) *     WWPJ(I,J) *(1.-WWPK(K))*WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)+1) + &
                  (1.-WWPI(I,J))*     WWPJ(I,J) *(1.-WWPK(K))*WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)+1) + &
                      WWPI(I,J) * (1.-WWPJ(I,J))*(1.-WWPK(K))*WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)+1) + &
                  (1.-WWPI(I,J))* (1.-WWPJ(I,J))*(1.-WWPK(K))*WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)+1) 

 SERR(I,J,K) =        WWPI(I,J) *     WWPJ(I,J) *    WWPK(K) *WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)  +NZLR_WW) + &
                  (1.-WWPI(I,J))*     WWPJ(I,J) *    WWPK(K) *WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)  +NZLR_WW) + &
                      WWPI(I,J) * (1.-WWPJ(I,J))*    WWPK(K) *WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)  +NZLR_WW) + &
                  (1.-WWPI(I,J))* (1.-WWPJ(I,J))*    WWPK(K) *WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)  +NZLR_WW) + &
                      WWPI(I,J) *     WWPJ(I,J) *(1.-WWPK(K))*WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)+1+NZLR_WW) + &
                  (1.-WWPI(I,J))*     WWPJ(I,J) *(1.-WWPK(K))*WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)+1+NZLR_WW) + &
                      WWPI(I,J) * (1.-WWPJ(I,J))*(1.-WWPK(K))*WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)+1+NZLR_WW) + &
                  (1.-WWPI(I,J))* (1.-WWPJ(I,J))*(1.-WWPK(K))*WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)+1+NZLR_WW) 

    ENDDO
  ENDDO
ENDDO

END SUBROUTINE SMPL_TL

SUBROUTINE SMPL_AD

IMPLICIT NONE

INTEGER(I4) :: I,J,K

DO K=GRD%KM,1,-1
  DO J=GRD%JM,1,-1
    DO I=GRD%IM,1,-1

    WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)  )=WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)  )+    WWPI(I,J) *WWPJ(I,J)     *    WWPK(K) *TERR_AD(I,J,K)
    WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)  )=WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)  )+(1.-WWPI(I,J))*WWPJ(I,J)     *    WWPK(K) *TERR_AD(I,J,K)
    WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)  )=WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)  )+    WWPI(I,J) *(1.-WWPJ(I,J))*    WWPK(K) *TERR_AD(I,J,K)
    WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)  )=WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)  )+(1.-WWPI(I,J))*(1.-WWPJ(I,J))*    WWPK(K) *TERR_AD(I,J,K)
    WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)+1)=WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)+1)+    WWPI(I,J) *WWPJ(I,J)     *(1.-WWPK(K))*TERR_AD(I,J,K)
    WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)+1)=WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)+1)+(1.-WWPI(I,J))*WWPJ(I,J)     *(1.-WWPK(K))*TERR_AD(I,J,K)
    WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)+1)=WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)+1)+    WWPI(I,J) *(1.-WWPJ(I,J))*(1.-WWPK(K))*TERR_AD(I,J,K)
    WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)+1)=WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)+1)+(1.-WWPI(I,J))*(1.-WWPJ(I,J))*(1.-WWPK(K))*TERR_AD(I,J,K)

    WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)  +NZLR_WW)=WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)  +NZLR_WW)+    WWPI(I,J) *WWPJ(I,J)     *    WWPK(K) *SERR_AD(I,J,K)
    WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)  +NZLR_WW)=WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)  +NZLR_WW)+(1.-WWPI(I,J))*WWPJ(I,J)     *    WWPK(K) *SERR_AD(I,J,K)
    WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)  +NZLR_WW)=WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)  +NZLR_WW)+    WWPI(I,J) *(1.-WWPJ(I,J))*    WWPK(K) *SERR_AD(I,J,K)
    WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)  +NZLR_WW)=WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)  +NZLR_WW)+(1.-WWPI(I,J))*(1.-WWPJ(I,J))*    WWPK(K) *SERR_AD(I,J,K)
    WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)+1+NZLR_WW)=WWT2(WWKI(I,J) ,WWKJ(I,J) ,WWKK(K)+1+NZLR_WW)+    WWPI(I,J) *WWPJ(I,J)     *(1.-WWPK(K))*SERR_AD(I,J,K)
    WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)+1+NZLR_WW)=WWT2(WWKI2(I,J),WWKJ(I,J) ,WWKK(K)+1+NZLR_WW)+(1.-WWPI(I,J))*WWPJ(I,J)     *(1.-WWPK(K))*SERR_AD(I,J,K)
    WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)+1+NZLR_WW)=WWT2(WWKI(I,J) ,WWKJ2(I,J),WWKK(K)+1+NZLR_WW)+    WWPI(I,J) *(1.-WWPJ(I,J))*(1.-WWPK(K))*SERR_AD(I,J,K)
    WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)+1+NZLR_WW)=WWT2(WWKI2(I,J),WWKJ2(I,J),WWKK(K)+1+NZLR_WW)+(1.-WWPI(I,J))*(1.-WWPJ(I,J))*(1.-WWPK(K))*SERR_AD(I,J,K)

    ENDDO
  ENDDO
ENDDO
END SUBROUTINE SMPL_AD

SUBROUTINE PRINT_MODERR(NIT)

IMPLICIT NONE
INTEGER(I4), INTENT(IN) :: NIT
INTEGER(I4) :: JI,JJ,JK,JT
CHARACTER(LEN=99) :: CF

WRITE(CF,'(A,I3.3,A)' ) 'wc4_iter_',NIT,'.nc'
CALL WRITE_VARNCDF(CF,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW,WW,'ww','degC/psu',&
'Model error')

WRITE(CF,'(A,I3.3,A)' ) 'wc4ad_iter_',NIT,'.nc'
CALL WRITE_VARNCDF(CF,NXLR_WW,NYLR_WW,2*NZLR_WW,NTS_WW,WW_AD,'ww_ad','degC/psu',&
'Model error')

END SUBROUTINE PRINT_MODERR

END MODULE WEAKLY
