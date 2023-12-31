SUBROUTINE ASSIGN_SLAERR

!.. OBSERVATIONAL ERRORS DEFINITION FOR SLA
!..
!.. ERRORS ARE DEFINED AS A SUM OF THE MDT ERROR
!.. PLUS
!.. A SATELLITE-DEPENDENT INSTRUMENTAL ERROR PLUS
!.. A LATLON-DEPENDENT REPRESENTIVENESS ERROR
!.. PLUS A LATITUDE-DEPENDENT OBSERVATION OPERATOR
!.. ERROR. A FACTOR TO ACCOUNT FOR TEMPORAL ERROR
!.. GROWTH AS A FUNCTION OF DISTANCE FROM ANALYSIS
!.. TIME IS NOT YET USED.
!..
!.. A.STORTO, 2009

USE SET_KND
USE GRD_STR
USE OBSDEF
USE OBS_STR
USE MYNETCDF
USE IOUNITS
USE RUN
USE MOD_GCDIST, ONLY : GCDIST
USE OBSHANDLING
USE MYFRTPROF
USE MYNETCDF
USE IEEE_ARITHMETIC
USE MPIREL

IMPLICIT NONE

REAL(R8) :: MDTERR,REPERR,INSERR,OOERR
REAL(R8) :: REP_ERR
REAL(R8) :: OO_ERR_MIN,OO_ERR_MAX,OO_ERR_LAT,IBERR
REAL(R8) :: MDT_ERRMAX,MDT_ERRMIN,MDT_ERRFIX
REAL(R8) :: GRD_MDTERR(GRD%IM,GRD%JM)
REAL(R8) :: ZMDTERR(GRD%IM,GRD%JM)
REAL(R8) :: ZREPERR(GRD%IM,GRD%JM)
REAL(R8) :: ZIBERR(GRD%IM,GRD%JM)
INTEGER(I4) :: NREPNAN(GRD%IM,GRD%JM)
REAL(R8), ALLOCATABLE :: ZLHAERR(:)
REAL(R8) :: REPERR_INFLAT
REAL(DP) :: TD, TD2, OETIM, ZT2
INTEGER(I4) :: K, STAT
INTEGER(I4) :: INTOT,ILATST,IINCR,JLAT
INTEGER(I4) :: KSATID, I, J, I2,KINTLAT,IAUX
LOGICAL :: LOUTPUT
REAL(R8) :: ZDIST
INTEGER(I4) :: NCID,IDVAR
INTEGER(I4) :: XIND1

#include "grids.h"

CALL MYFRTPROF_WALL('ASSIGN_SLAERR: ASSIGN SLA OBS ERROR',0)

!===================================================================
!    DEFINE INSTRUMENTAL AND REPRESENTIVENESS ERRORS CONSTANTS

REP_ERR=0.03_R8

OO_ERR_MIN=0.01_R8
OO_ERR_MAX=0.05_R8
OO_ERR_LAT=17.5_R8

MDT_ERRMAX = 0.1_R8
MDT_ERRMIN = 0.001_R8
MDT_ERRFIX = 0.02_R8

REPERR_INFLAT = SLA_REPERR_INFLAT

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) ' *** INFLATION FOR SLA REPRESENTATIVENESS ERROR ***'
WRITE(IOUNLOG,*) ' INFLATION COEFFICIENT: ', REPERR_INFLAT
WRITE(IOUNLOG,*)
CALL FLUSH(IOUNLOG)

!
!===================================================================

LOUTPUT=.FALSE.

!... ASSIGN ERROR FOR SLA OBSERVATIONS

! MDT error
CALL GETNCVAR('MDT.nc','mdt_err',GRD%IM,GRD%JM,ZMDTERR)
ZMDTERR = ABS(ZMDTERR)

IF(LLVARMDT .OR. NCONF .EQ. 202) THEN
   GRD%CMDT_STDEV = ABS( ZMDTERR )
   IF( LLVARMDT ) WHERE ( GRD%CMDT_STDEV .GT. 0.02_R8 ) GRD%CMDT_STDEV = 0.02_R8
   GRD%CMDT_STDEV = GRD%CMDT_STDEV / ZVARMDT_RATIOSD
ENDIF

! SLA representativity error
CALL GETNCVAR('SLA_repres.nc','stde',GRD%IM,GRD%JM,ZREPERR)

! Inverse barometer error
IF( LL_IBERR ) THEN
  WRITE(IOUNLOG,*) ' READING IB ERROR FOR SLA'
  CALL GETNCVAR('IB_error.nc','iberror',GRD%IM,GRD%JM,ZIBERR)
ELSE
  WRITE(IOUNLOG,*) ' IB ERROR FOR SLA IS NOT ACCOUNTED FOR'
ENDIF

! Observation operator error (function of latitude)
IF(LHA_ERR) THEN
  OPEN(751,FILE='lha_err.dat',STATUS='OLD',IOSTAT=STAT)
  IF(STAT.NE.0) CALL ABOR1('ASSIGN_SLAERR: CANNOT OPEN LHA_ERR.DAT')
  READ(751,*,IOSTAT=STAT) ILATST,IINCR,INTOT
  IF(STAT.NE.0) CALL ABOR1('ASSIGN_SLAERR: CANNOT READ LHA_ERR.DAT')
  ALLOCATE(ZLHAERR(INTOT))
  DO JLAT=1,INTOT
     READ(751,*,IOSTAT=STAT) IAUX, ZLHAERR(JLAT)
     IF(STAT.NE.0) CALL ABOR1('ASSIGN_SLAERR: CANNOT READ LHA_ERR.DAT')
  ENDDO
  CLOSE(751)
ENDIF

IF(LOUTPUT) THEN
    OPEN (91,FILE='slaerr.dat.'//TRIM(CMPIDOM))
    WRITE(91,*) '  OBSN     TRACK      LAT MDTER REPER INSER OBOER TOTER'
ENDIF

CYOBS : DO K=1,SLA%NO

   IF(SLA%FLC(K).EQ.0 ) CYCLE CYOBS

   KSATID=SLA%KSAT(K)

   IF(KSATID.EQ.0) THEN
        WRITE(IOUNERR,*) 'TRACK/SAT IS ', SLA%TRACK(K),KSATID
        CALL ABOR1('UNKNOWN SATELLITE TRACK IN ASSIGN_SLAERR')
   ENDIF

   IF(LLUSEMDTERR) THEN
     MDTERR= OSUM(SLA%PQ(K,:) , ZMDTERR(SLA%IB(K,:),SLA%JB(K,:)) ) 
     MDTERR = MIN(MAX(MDTERR,0.008_R8),0.5_R8)
   ELSE
     MDTERR=0._R8
   ENDIF

   IF(LL_IBERR) THEN
     IBERR= OSUM(SLA%PQ(K,:) , ZIBERR(SLA%IB(K,:),SLA%JB(K,:)) ) 
     IBERR = MIN(MAX(IBERR,0._R8),0.05_R8)
   ELSE
     IBERR=0._R8
   ENDIF

   !... REPRESENTATIVENESS ERROR
   REPERR= OSUM(SLA%PQ(K,:) , ZREPERR(SLA%IB(K,:),SLA%JB(K,:)) ) 

   REPERR=REPERR_INFLAT*MAX(MIN(REPERR,0.5_R8),0.005_R8)

   !... INSTRUMENTAL ERROR
   INSERR=SLA_ERR(KSATID)

   !... OBS OPERATOR (LHA) ERROR
   IF(LHA_ERR) THEN
    KINTLAT=NINT( (SLA%LAT(K)-REAL(ILATST))/REAL(IINCR) )+1
    IF(KINTLAT.LT.1.OR.KINTLAT.GT.INTOT) THEN
      WRITE(IOUNERR,*) 'LATITUDE=',SLA%LAT(K),',  KINTLAT=',KINTLAT
      WRITE(IOUNERR,*) 'UNSUPPORTED LHA ERROR'
      CALL ABOR1('ASSIGN_SLA_ERR: LHA ERROR UNDETERMINED')
    ENDIF
    OOERR=ZLHAERR(KINTLAT)
   ELSE
    OOERR=0._R8
   ENDIF

!... SUM ERRORS
   !SLA%ERR(K) = SQRT(MDTERR*MDTERR + REPERR*REPERR + INSERR*INSERR + OOERR*OOERR + IBERR*IBERR)
   SLA%ERR(K) = INSERR

!... INFLAT FOR TIME DISTANCE
   TD=MIN(ZTIME_WEIGTH,ABS(SLA%TDIST(K)))
   TD2=TD*TD
   ZT2=ZTIME_WEIGTH*ZTIME_WEIGTH
   OETIM = 1._DP / EXP( -TD2 / ZT2 )
   SLA%ERR(K) = SLA%ERR(K)*REAL(OETIM,R8)
   
   IF(LOUTPUT) THEN
    WRITE(91,'(I6,X,I6,X,F8.3,5(X,F6.3))') K, SLA%TRACK(K),&
    & SLA%LAT(K),MDTERR,REPERR,INSERR,OOERR,SLA%ERR(K)
   ENDIF

ENDDO CYOBS

IF(LOUTPUT) CLOSE(91)

CALL MYFRTPROF_WALL('ASSIGN_SLAERR: ASSIGN SLA OBS ERROR',1)

END SUBROUTINE ASSIGN_SLAERR
