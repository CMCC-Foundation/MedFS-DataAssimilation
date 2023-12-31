SUBROUTINE SFILTATSLAL(NO,ZLON,ZLAT,ZTIM,ZVAL,NKN,ZCL)

!.. PROCESS ALONG-TRACK SLA TO FILTER OUT HIGH FREQUENCIES
!.. THROUGH A LOW-PASS LANCZOS FILTER
!.. A.S.

USE SET_KND
USE FILTRW
USE IOUNITS
USE OCEANTOOLS, ONLY : DIST
USE MOD_GCDIST

IMPLICIT NONE

INTEGER(I4) , INTENT(IN) :: NO,NKN
REAL(R8), DIMENSION(NO),INTENT(IN) :: ZLON,ZLAT,ZTIM
REAL(R8), INTENT(IN) :: ZCL
REAL(R8), DIMENSION(NO),INTENT(INOUT) :: ZVAL
REAL(R8), DIMENSION(NO) :: ZWRK
INTEGER(I4) :: JI,JJ,JS,JE,NOFF
REAL(R8), ALLOCATABLE :: ZCOE(:)
INTEGER(I4) :: KC,IERR,MAXC,ILUT
REAL(R8) :: ZDIST,ZAUX,ZCOF,ZTOT,ZDD
REAL(R8),ALLOCATABLE :: WEI1D(:),WEI2D(:,:),RES(:),FRE(:)
REAL(R8), PARAMETER :: ZSLARES = 20._R8
LOGICAL :: LLLANCVAR,LLPP
INTEGER(I4) :: JJR,ICO,NRESLAT
INTEGER(I4),PARAMETER :: MAXRES=200
REAL(R8) :: ZRESS(MAXRES),ZRESM(MAXRES),ZFREQ(MAXRES),ZTABLAT(MAXRES)

ZWRK=ZVAL

LLLANCVAR=.TRUE.
LLPP=.FALSE.
NRESLAT=81

IF(LLLANCVAR) THEN
   OPEN(91,FILE='FILRES.DAT',STATUS='OLD')
   DO JJR=1,NRESLAT
      READ(91,*) ZTABLAT(JJR),ZRESS(JJR),ZRESM(JJR),ZFREQ(JJR)
   ENDDO
   CLOSE(91)
   ALLOCATE(WEI2D(NRESLAT,NKN))
   DO JJR=1,NRESLAT
      ZCOF=ZFREQ(JJR)
      CALL DFILTRQ(NKN,ZCOF,0._R8,1._R8,0,WEI2D(JJR,:),RES,FRE,IERR)
   ENDDO
ELSE
   ALLOCATE(WEI1D(NKN),RES(2*NKN-1),FRE(2*NKN-1))
   ZCOF=ZSLARES/ZCL
   CALL DFILTRQ(NKN,ZCOF,0._R8,1._R8,0,WEI1D,RES,FRE,IERR)
   IF(IERR.NE.0) THEN
      WRITE(IOUNERR,*) ' DFILTRQ IN SFILTATSLAL RETURNED',IERR
      CALL ABOR1('FILTER WEIGHTS COMPUTATION FAILED')
   ENDIF
ENDIF

MAXC=NKN
ALLOCATE(ZCOE(MAXC))

NOFF=(NKN-1)/2
ZVAL=0._R8

!.. SMOOTHED BOUNDARY CONDITIONS ARE APPLIED

WRITE(1027,*) '///'

CYO : DO JI=1,NO

  IF(ABS(ZWRK(JI)).GT.1000._R8) THEN
    ZVAL(JI)=ZWRK(JI)
    CYCLE CYO
  ENDIF

  ZTOT=0._R8
  ICO=0
  JS=MAX(1,JI-NOFF)
  JE=MIN(NO,JI+NOFF)

  KC=0
  ILUT=NINT( (ZLAT(JI) + 80._R8)/2._R8 + 1._R8 )
  DO JJ=JS,JE
    KC=KC+1
    IF(ABS(ZWRK(JJ)).LT.1000._R8) THEN
       CALL GCDIST(ZLAT(JI),ZLON(JI),ZLAT(JJ),ZLON(JJ),ZDD)
       IF(ZDD.LT.ZRESS(ILUT)*(REAL(ABS(JJ-JI),KIND=R8)+0.5_R8) ) THEN
         ZVAL(JI)=ZVAL(JI)+WEI2D(ILUT,KC)*ZWRK(JJ)
         ZTOT=ZTOT+WEI2D(ILUT,KC)
         ICO=ICO+1
       ENDIF
    ENDIF
  ENDDO
  IF(ICO.EQ.0) ZVAL(JI)=ZWRK(JI)

  ZVAL(JI)=ZVAL(JI)/ZTOT

  IF(JI.LT.NO.AND.LLPP) THEN
    CALL GCDIST(ZLAT(JI),ZLON(JI),ZLAT(JI+1),ZLON(JI+1),ZDD)
    WRITE(1027,'(2I5,3F9.2,X,2F12.8,I3)') &
    & JI,NO,ZLON(JI),ZLAT(JI),ZDD,ZWRK(JI),ZVAL(JI),ICO
  ENDIF

ENDDO CYO

RETURN
END SUBROUTINE SFILTATSLAL
