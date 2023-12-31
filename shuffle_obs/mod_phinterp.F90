! Define by CPP macro the dist function (default to DIST2)
#if defined DIST_FUNCTION
#define _MYDIST_ DIST_FUNCTION
#else
#define _MYDIST_ HAVERSINE
#endif
MODULE PHINTERP

!
! INTERPOLATION PARAMETER
!

USE SET_KND, ONLY : R8, I4
USE MYNETCDF

IMPLICIT NONE

INTEGER, PARAMETER :: IOUNOUT=6
INTEGER, PARAMETER :: IOUNLOG=IOUNOUT

TYPE GRID_T
        INTEGER(I4)              ::  IM, JM
        REAL(R8),ALLOCATABLE     :: LON(:,:), LAT(:,:)
END TYPE GRID_T

TYPE (GRID_T)                 :: GRD

REAL(R8), PARAMETER :: ZLATB = 20._R8, EPSLN=1.E-12_R8
REAL(R8) :: DEG2RAD, ZPI, MAXLAT
INTEGER(I4), PARAMETER :: NNSEC = 10, NITER=7
INTEGER(I4), PARAMETER :: ISTP = 2**(NITER-1)
INTEGER(I4) :: ICM
INTEGER(I4), PARAMETER :: NLS = 180._R8 / ZLATB
REAL(R8), ALLOCATABLE, PRIVATE :: INTLON(:,:)
INTEGER(I4), PRIVATE :: JAST(NLS), JAEN(NLS)
REAL(R8) :: MDIST(NLS,2), ZERAD(NLS)
LOGICAL :: LLINIT = .FALSE.

CONTAINS

SUBROUTINE PREPINTERP_INIT

!
! BUILD-UP TABLE FOR SPEEDING UP PREPINTERP
!

IMPLICIT NONE

INTEGER(I4) :: I, KTMP(1), JI, JJ, ISTP, J1, J2
REAL(R8)    :: LAT0, LAT1, MINLA(GRD%JM), MAXLA(GRD%JM), LAT0B, LAT1B
REAL(R8),PARAMETER   :: ZA1 = 6378137.0_R8
REAL(R8),PARAMETER   :: ZB1 = 6356752.3_R8

WRITE(IOUNLOG,*) 'PINT_START' ; CALL FLUSH(IOUNLOG)

CALL GETNCDIM('GRID_CR.nc','x',GRD%IM)
CALL GETNCDIM('GRID_CR.nc','y',GRD%JM)

WRITE(IOUNOUT,*) ' DIMS ARE ',GRD%IM, GRD%JM

ALLOCATE( GRD%LON(GRD%IM, GRD%JM) )
ALLOCATE( GRD%LAT(GRD%IM, GRD%JM) )

CALL GETNCVAR ('GRID_CR.nc','lon',GRD%IM,GRD%JM,GRD%LON)
CALL GETNCVAR ('GRID_CR.nc','lat',GRD%IM,GRD%JM,GRD%LAT)

ZPI = 2._R8 * ASIN(1._R8)
DEG2RAD = ZPI/180._R8

ICM = GRD%IM-1 + 2*ISTP
ALLOCATE( INTLON(ICM,GRD%JM) )

DO JJ=1,GRD%JM
  DO JI=1,GRD%IM
     INTLON(JI,JJ) = GRD%LON(JI,JJ)
  ENDDO
  DO JI=1,2*ISTP
     INTLON(GRD%IM-1+JI,JJ) = GRD%LON(JI+1,JJ)
  ENDDO
ENDDO

MAXLAT= MAXVAL(GRD%LAT)

DO I=1,NLS
   LAT0=-90._R8+(I-1)*ZLATB
   LAT1=-90._R8+(I)*ZLATB

   LAT0B=MAX(-90._R8,LAT0-2._R8)
   LAT1B=MIN( 90._R8,LAT1+2._R8)
   J1=GRD%JM
   J2=1
   DO JJ=1,GRD%JM
      IF( ANY(GRD%LAT(:,JJ) .GE. LAT0B .AND. GRD%LAT(:,JJ) .LE. LAT1B ) ) THEN
          J1=MIN(J1,JJ)
          J2=MAX(J2,JJ)
      ENDIF 
   ENDDO
   J1 = MAX(1,J1-NNSEC)
   J2 = MIN(GRD%JM,J2+NNSEC)
   JAST(I) = MIN(J1,J2)
   JAEN(I) = MAX(J1,J2)

   ZERAD(I) = 0.5_R8*SQRT( ZA1**4*COS(LAT0)**2 + ZB1**4*SIN(LAT0)**2 ) / &
                   & SQRT( ZA1**2*COS(LAT0)**2 + ZB1**2*SIN(LAT0)**2 ) + &
              0.5_R8*SQRT( ZA1**4*COS(LAT1)**2 + ZB1**4*SIN(LAT1)**2 ) / &
                   & SQRT( ZA1**2*COS(LAT1)**2 + ZB1**2*SIN(LAT1)**2 )
   MDIST(I,1) = MAXVAL( (GRD%LON(2:GRD%IM,JAST(I):JAEN(I))-&
                         GRD%LON(1:GRD%IM-1,JAST(I):JAEN(I))) )!! &
   MDIST(I,1) = 2.5_R8
   MDIST(I,2) = MAXVAL( (GRD%LAT(:,JAST(I)+1:JAEN(I))-&
                         GRD%LAT(:,JAST(I):JAEN(I)-1)) )
   WRITE(*,*) ' LATBOX ', I, ' LIMITS : ',LAT0,LAT1,' MIN/MAX ', JAST(I), JAEN(I), &
   & MINLA(JAST(I)),MAXLA(JAST(I)), MINLA(JAEN(I)),MAXLA(JAEN(I))
   WRITE(*,*) '       MAXD ',MDIST(I,1:2), ' RADIUS :',ZERAD(I)
ENDDO

WRITE(IOUNLOG,*) 'PINT_END' ; CALL FLUSH(IOUNLOG)

END SUBROUTINE PREPINTERP_INIT

SUBROUTINE PREPAUX(LO,LA,IS,IE,JS,JE,MDS,ZR)

  IMPLICIT NONE

  REAL(R8), INTENT(IN) :: LO, LA
  INTEGER(I4), INTENT(OUT):: IS, IE, JS, JE
  REAL(R8), INTENT(OUT):: MDS(2),ZR

  INTEGER(I4) :: KJ

  KJ=INT((LA+90._R8)/ZLATB)+1
  JS=JAST(KJ)
  JE=JAEN(KJ)
  MDS(1:2)=MDIST(KJ,1:2)
  IS=1
  IE=ICM
  ZR=ZERAD(KJ)

END SUBROUTINE PREPAUX

SUBROUTINE PREPINTERP(LON,LAT,IBS,JBS,PBS)

IMPLICIT NONE

REAL(R8)   , INTENT(IN)  :: LON, LAT
INTEGER(I4), INTENT(OUT) :: IBS(4), JBS(4)
REAL(R8)   , INTENT(OUT) :: PBS(4)

INTEGER(I4), PARAMETER   :: NPP = 25
INTEGER(I4) :: JI, JJ, K, KK, I0, J0, ITM(NPP), JTM(NPP), ISN(NPP), KI
REAL(R8)    :: ZMD(2), ZDIST, ZD, ZDA(NPP), ZRD
INTEGER(I4) :: ISTP, IST, IEN, JST, JEN, ISTART, IEND, JSTART, JEND, ITER, IOS(5)

IF( .NOT.LLINIT ) THEN
   CALL PREPINTERP_INIT
   LLINIT = .TRUE.
ENDIF

CALL PREPAUX(LON,LAT,ISTART,IEND,JSTART,JEND,ZMD,ZRD)

ZDIST=1.E+20_R8
ISTP = 2**(NITER-1)
IST = ISTART
JST = JSTART
IEN = IEND
JEN = JEND
I0=0
J0=0

IF(LAT.GT.MAXLAT) RETURN

ISTP=1
!!!ISTP=2
!DO ITER=1,2
  DO JJ=JST,JEN,ISTP
    DO JI=IST,IEN,ISTP

      ZD = _MYDIST_(LON,LAT,INTLON(JI,JJ),GRD%LAT(JI,JJ))*ZRD
      IF (ZD .LT. ZDIST) THEN
         I0=JI
	 J0=JJ
	 ZDIST=ZD
      ENDIF

    ENDDO
  ENDDO
!  ISTP=1
!  IST=MAX(I0-9,1)
!  IEN=MIN(ICM,I0+9)
!  JST=MAX(JSTART,J0+9)
!  JEN=MIN(JEND,J0+9)
!ENDDO

IF(I0.EQ.0 .OR. J0.EQ.0) RETURN

  IF( HELPLD(LON,INTLON(I0,J0)) .GT. ZMD(1) .OR. &
      ABS(LAT-GRD%LAT(I0,J0)) .GT. ZMD(2) ) THEN
            I0 = 0
            J0 = 0
            RETURN
  ENDIF

IF(I0.GT.GRD%IM) I0=I0-GRD%IM+2

IOS=(/I0-2,I0-1,I0,I0+1,I0+2/)
WHERE( IOS .LT. 1      ) IOS = IOS+GRD%IM-2
WHERE( IOS .GT. GRD%IM ) IOS = IOS-GRD%IM+2

K=0
DO JJ=MAX(J0-2,1),MIN(GRD%JM,J0+2)
  DO KI=1,5
     JI=IOS(KI)
     K=K+1
     ZDA(K)=_MYDIST_(LON,LAT,INTLON(JI,JJ),GRD%LAT(JI,JJ))
     ITM(K)=JI
     JTM(K)=JJ
  ENDDO
ENDDO

CALL INDEX_SORT(ZDA(1:K),ISN(1:K),K)

DO K=1,4
   KK = ISN(K)
   IBS(K) = ITM(KK)
   JBS(K) = JTM(KK)
   PBS(K) = ZDA(KK)
ENDDO
PBS = PBS / SUM(PBS)

!* WE PASS DATA FROM EAST TO WEST
WHERE(IBS.EQ.1 .OR. IBS.EQ.2) IBS = IBS+GRD%IM-2

END SUBROUTINE PREPINTERP

REAL(R8) FUNCTION DIST(LO1,LA1,LO2,LA2)
IMPLICIT NONE

REAL(R8), INTENT(IN) :: LO1,LA1,LO2,LA2

DIST = ( (LO1-LO2)*cos( (LA1+LA2)*zpi/180.*0.5 ) )**2 + (LA1-LA2)**2

END FUNCTION DIST

REAL(R8) FUNCTION DIST2(LO1,LA1,LO2,LA2)
IMPLICIT NONE

REAL(R8), INTENT(IN) :: LO1,LA1,LO2,LA2
REAL(R8) :: R1(3), R2(3), CROSS(3), s, dot, ang, LOA2, LOA1

LOA1=LO1+EPSLN
LOA2=LO2+EPSLN
R1(1) = COS(DEG2RAD*LOA1)*COS(DEG2RAD*LA1)
R1(2) = SIN(DEG2RAD*LOA1)*COS(DEG2RAD*LA1)
R1(3) = SIN(DEG2RAD*LA1)
R2(1) = COS(DEG2RAD*LOA2)*COS(DEG2RAD*LA2)
R2(2) = SIN(DEG2RAD*LOA2)*COS(DEG2RAD*LA2)
R2(3) = SIN(DEG2RAD*LA2)

cross(1) = r1(2)*r2(3)-r1(3)*r2(2)
cross(2) = r1(3)*r2(1)-r1(1)*r2(3)
cross(3) = r1(1)*r2(2)-r1(2)*r2(1)

s = sqrt(cross(1)**2+cross(2)**2+cross(3)**2)

s = min(s,1._R8-epsln)

dot = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)

if (dot > 0) then
    ang = asin(s)
else if (dot < 0) then
    ang = zpi - asin(s)
else
    ang = zpi/2._R8
endif

DIST2 = ABS(ANG)

END FUNCTION DIST2

REAL(R8) FUNCTION GCDIST(LO1,LA1,LO2,LA2)

IMPLICIT NONE

REAL(R8),INTENT(IN)  :: LO1,LA1,LO2,LA2
REAL(R8) ::  RPD,PI,PI2T,DPR,PIO2
REAL(R8)             :: ZLATS,ZLONS,ZLATE,ZLONE
REAL(R8)             :: ZW,ZC1,ZC2,LL,DLO

ZLATS=LA1
ZLATE=LA2
ZLONS=LO1
ZLONE=LO2

PI2T = 2._R8*ZPI
PI2T = ZPI*0.5_R8
RPD=ZPI/180._R8
DPR=180._R8/ZPI

IF(ZLONS.LT.0._R8) ZLONS = ZLONS+360._R8
IF(ZLONE.LT.0._R8) ZLONE = ZLONE+360._R8

DLO = (ZLONE-ZLONS)*RPD

IF(DLO.LT.-PI) DLO=DLO+PI2T
IF(DLO.GT.PI)  DLO=DLO-PI2T

ZLATS=ZLATS*RPD
ZLATE=ZLATE*RPD

IF(ZLATS.EQ.ZLATE .AND. DLO .EQ. 0._R8) THEN
   GCDIST=0._R8
   RETURN
ENDIF

LL = ZLATE - ZLATS

ZC1 = PIO2 - ZLATS
ZC2 = PIO2 - ZLATE

GCDIST = ((1._R8-COS(DLO))*0.5_R8)*COS(ZLATS)*COS(ZLATE) +&
     & ((1._R8-COS(LL))*0.5_R8)

ZW = 1._R8 - 2._R8*GCDIST
IF( ZW .GT. 1._R8 ) ZW = 1._R8
IF( ZW .LT.-1._R8 ) ZW =-1._R8

GCDIST = ACOS(ZW)*DPR*111.12000071117_R8

END FUNCTION GCDIST

REAL(R8) FUNCTION HELPLD (LO1,LO2)
IMPLICIT NONE

REAL(R8),INTENT(IN)  :: LO1,LO2
REAL(R8) :: DLO

DLO=ABS(LO1-LO2)
IF( DLO .GT. 180._R8) THEN
    IF( LO1.GT.0._R8 .AND.LO2.LT.0 ) THEN
        DLO = ABS( LO1-(LO2+360._R8) )
    ELSE
        DLO = ABS( LO1+360._R8-LO2 )
    ENDIF
ENDIF

HELPLD=DLO

END FUNCTION HELPLD

REAL(R8) FUNCTION HAVERSINE(LO1,LA1,LO2,LA2)

! A.S. 20.01.2010
! Haversine formula to compute distance (km)
! between points in latlon coordinates
!
! Use latitudinally-dependent Earth radius

IMPLICIT NONE

REAL(R8),INTENT(IN)  :: LO1,LA1,LO2,LA2
REAL(R8)             :: ZLATS,ZLONS,ZLATE,ZLONE,DLAT,DLON
REAL(R8)             :: ZA, ZC, ZR
REAL(R8),PARAMETER   :: ZA1 = 6378.1370_8
REAL(R8),PARAMETER   :: ZB1 = 6356.7523_8

ZLATS=LA1*DEG2RAD
ZLATE=LA2*DEG2RAD
ZLONS=LO1*DEG2RAD
ZLONE=LO2*DEG2RAD

DLAT=ABS(ZLATS-ZLATE)
DLON=ABS(ZLONS-ZLONE)

ZA = SIN( DLAT/2._8 )*SIN( DLAT/2._8 ) + COS(ZLATS)*COS(ZLATE)*&
& SIN(DLON/2._8)*SIN(DLON/2._8)
ZC = 2._8 * ATAN2(SQRT(ZA),SQRT(1._8-ZA))
ZR = SQRT( ZA1**4*COS(ZLATS)**2 + ZB1**4*SIN(ZLATS)**2 ) / &
   & SQRT( ZA1**2*COS(ZLATS)**2 + ZB1**2*SIN(ZLATS)**2 )
HAVERSINE = ZC*ZR

END FUNCTION HAVERSINE
END MODULE PHINTERP
