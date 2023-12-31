SUBROUTINE ASSIGN_ICOERR

!.. OBS ERRORS ASSIGNEMENT FOR IN-SITU OBSERVATIONS
!..
!.. ANDREA STORTO, 2009


USE SET_KND
USE OBS_STR
USE OBSDEF
USE NETCDF
USE OBSHANDLING, ONLY : ZTIME_WEIGTH, ICOSST_DEF_ERR
USE MYFRTPROF
USE MPIREL, ONLY : MYPROC

IMPLICIT NONE

REAL (R8) :: OERR,RDEP,OEREP
REAL (DP) :: OETIM
INTEGER(I4) :: JK,ITYP,IPAR,IDEP
INTEGER(I4) :: NDEPTHS
INTEGER(I4) :: JJX,JJY
INTEGER(I4) :: NCID, IDVAR, IDDIM
REAL (R8) :: ZDEP, ZLON, ZLAT
REAL (DP) :: TD, TD2, ZT2
REAL (R8), ALLOCATABLE :: ICOERRS(:,:), INSDEP(:)
REAL (R8) :: INSFAC(4)
INTEGER(I4), PARAMETER :: NXR=1440, NYR=720
REAL (R8) :: REPERR(NXR,NYR)

REAL(R8), PARAMETER :: MIN_ERR = 0.4_R8

#include "obs_events.h"

 CALL MYFRTPROF_WALL('ASSIGN_ICOERR: ASSIGN ICOADS OBS ERROR',0)

!... READ FROM FILE

CALL CHECK( NF90_OPEN('insitu_errors.nc',NF90_NOWRITE,NCID) )
CALL CHECK( NF90_INQ_VARID (NCID, 'reperr', IDVAR) )
CALL CHECK( NF90_GET_VAR (NCID,IDVAR,REPERR) )
CALL CHECK( NF90_CLOSE (NCID) )

!... FOR INS, CYCLE THE OBS

CYOBS : DO JK=1,INS%NO

     IF(INS%FLC(JK) .NE. 1 ) CYCLE CYOBS

       ! REVISION OF OBS ERRORS

       IF( INS%ERR(JK) .LT. 0._R8   ) INS%ERR(JK) = ICOSST_DEF_ERR
       IF( INS%ERR(JK) .LT. MIN_ERR ) INS%ERR(JK) = MIN_ERR

       TD=MIN(ZTIME_WEIGTH,ABS(INS%TDIST(JK)))
       TD2=TD*TD
       ZT2=ZTIME_WEIGTH*ZTIME_WEIGTH

       ZLON=INS%LON(JK)
       ZLAT=INS%LAT(JK)
       IF( ZLON.LT.0._R8) ZLON=ZLON+360._R8
       JJX=INT( ZLON/0.25_R8 )+1
       JJY=INT( (ZLAT+90._R8)/0.25_R8 )+1

       IF( JJX.GT.NXR .OR. JJX.LT.1 .OR. &
           JJY.GT.NYR .OR. JJY.LT.1 ) THEN
         INS%FLC(JK) = 0
         INS%EVE(JK) = KEVE_AOER
       ELSE
         OETIM = 1._DP / EXP( -TD2 / ZT2 )
         OEREP=MIN(MAX(REPERR(JJX,JJY),1._R8),2._R8)
         OERR = INS%ERR(JK)
         INS%ERR(JK) = OERR*REAL(OETIM,R8)*OEREP
       ENDIF

ENDDO CYOBS

 CALL MYFRTPROF_WALL('ASSIGN_ICOERR: ASSIGN ICOADS OBS ERROR',1)
END SUBROUTINE ASSIGN_ICOERR
