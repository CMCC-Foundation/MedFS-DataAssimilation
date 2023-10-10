SUBROUTINE OBS_VEC

!-----------------------------------------------------------------------
!                                                                      !
! CREATE THE OBSERVATIONAL VECTOR                                      !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE OBS_STR
 USE IOUNITS, ONLY : IOUNERR, IOUNLOG, IOUNOUT
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL

 IMPLICIT NONE

  INTEGER(I4)    ::  K, I

CALL MYFRTPROF_WALL('OBS_VEC: FORM OBS VECTOR',0)

! -------
! DEFINE OBSERVATIONAL VECTOR

    OBS%NO = SLA%NC + ARG%NC + XBT%NC + GLD%NC + VEL%NC + TRJ%NC

    WRITE(IOUNOUT,*) ' NUMBER OF OBSERVATIONS: ', OBS%NO
    WRITE(IOUNLOG,*) ' NUMBER OF OBSERVATIONS: ', OBS%NO

   IF (OBS%NO .LT. 1 ) THEN
      CALL ABOR1('NO OBSERVATIONS')
   ENDIF

   ALLOCATE ( OBS%INC(OBS%NO), OBS%AMO(OBS%NO), OBS%RES(OBS%NO))
   ALLOCATE ( OBS%ERR(OBS%NO), OBS%GRA(OBS%NO))

   K=0

! SLA OBSERVATIONS
   DO I=1,SLA%NO
    IF(SLA%FLC(I).EQ.1)THEN
     K=K+1
     OBS%RES(K) = SLA%RES(I)
     OBS%ERR(K) = SLA%ERR(I)
    ENDIF
   ENDDO

! ARGO OBSERVATIONS
   DO I=1,ARG%NO
    IF(ARG%FLC(I).EQ.1)THEN
     K=K+1
     OBS%RES(K) = ARG%RES(I)
     OBS%ERR(K) = ARG%ERR(I)
    ENDIF
   ENDDO

! XBT OBSERVATIONS
   DO I=1,XBT%NO
    IF(XBT%FLC(I).EQ.1)THEN
     K=K+1
     OBS%RES(K) = XBT%RES(I)
     OBS%ERR(K) = XBT%ERR(I)
    ENDIF
   ENDDO

! GLIDER OBSERVATIONS
   DO I=1,GLD%NO
    IF(GLD%FLC(I).EQ.1)THEN
     K=K+1
     OBS%RES(K) = GLD%RES(I)
     OBS%ERR(K) = GLD%ERR(I)
    ENDIF
   ENDDO

! VELOCITY OBSERVATIONS
   DO I=1,VEL%NO
    IF(VEL%FLC(I).EQ.1)THEN
     K=K+1
     OBS%RES(K) = VEL%RES(I)
     OBS%ERR(K) = VEL%ERR(I)
    ENDIF
   ENDDO

! VELOCITY OBSERVATIONS
   DO I=1,TRJ%NO
    IF(TRJ%FLC(I).EQ.1)THEN
     K=K+1
     OBS%RES(K) = TRJ%RES(I)
     OBS%ERR(K) = TRJ%ERR(I)
    ENDIF
   ENDDO

CALL MYFRTPROF_WALL('OBS_VEC: FORM OBS VECTOR',1)
END SUBROUTINE OBS_VEC
