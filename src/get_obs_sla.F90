SUBROUTINE GET_OBS_SLA
#ifdef MFS
!-----------------------------------------------------------------------
!                                                                      !
! LOAD SLA OBSERVATIONS                                                !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

  INTEGER(I4)   ::  K, N_TOPEX, N_ERS
  INTEGER(I4)   ::  I1, J1, KK, I, J, ITER
  REAL(R8)      ::  P1, Q1, P, Q, SUMT, SUMI
  REAL(DP)      ::  TIMP

  SLA%NO = 0
  SLA%NC = 0

  IF(OBS%SLA.EQ.0)RETURN

  OPEN(511,FILE='SLA_MIS.DAT',FORM='UNFORMATTED',STATUS='OLD')

! ---
! ALLOCATE MEMORY FOR OBSERVATIONS

   READ(511) SLA%NO
   PRINT*,'NUMBER OF SLA OBSERVATIONS: ',  SLA%NO

   IF(SLA%NO.EQ.0)RETURN

   ALLOCATE ( SLA%INO(SLA%NO), SLA%FLG(SLA%NO), SLA%FLC(SLA%NO))
   ALLOCATE ( SLA%LON(SLA%NO), SLA%LAT(SLA%NO), SLA%TIM(SLA%NO))
   ALLOCATE ( SLA%VAL(SLA%NO), SLA%BAC(SLA%NO), SLA%INC(SLA%NO))
   ALLOCATE ( SLA%BIA(SLA%NO), SLA%ERR(SLA%NO))
   ALLOCATE ( SLA%RES(SLA%NO), SLA%B_A(SLA%NO))
   ALLOCATE ( SLA%IB(SLA%NO), SLA%JB(SLA%NO))
   ALLOCATE ( SLA%PB(SLA%NO), SLA%QB(SLA%NO))
   ALLOCATE ( SLA%PQ1(SLA%NO), SLA%PQ2(SLA%NO), SLA%PQ3(SLA%NO), SLA%PQ4(SLA%NO))
   ALLOCATE ( SLA%DPT(SLA%NO))

! ---
! INITIALISE QUALITY FLAG
   SLA%FLG(:) = 1
   SLA%FLC(:) = 1

! ---
! LEVEL CORRESPONDING TO THE MINIMUM DEPTH
   SLA%KDP=GRD%KM
  DO K=GRD%KM, 1, -1
   IF(GRD%DEP(K).GE.SLA%DEP) SLA%KDP = K
  ENDDO

   PRINT*,' SLA%DEP: ',SLA%KDP

       READ (511)                                              &
        SLA%INO(1:SLA%NO), SLA%FLG(1:SLA%NO)                   &
       ,SLA%LON(1:SLA%NO), SLA%LAT(1:SLA%NO), SLA%TIM(1:SLA%NO)&
       ,SLA%VAL(1:SLA%NO), SLA%BAC(1:SLA%NO)                   &
       ,SLA%ERR(1:SLA%NO), SLA%RES(1:SLA%NO)                   &
       ,SLA%IB(1:SLA%NO), SLA%JB(1:SLA%NO)                     &
       ,SLA%PB(1:SLA%NO), SLA%QB(1:SLA%NO)
    CLOSE(511)

! ---
! REMOVE BIAS ALONG EACH TRACK AND OBSERAVTAIONS WITH LARGE RESIDUALS

 DO ITER = 1,3

!BIAS
   TIMP = SLA%TIM(1)
   I1 = 1
 DO K=2,SLA%NO
  IF(SLA%TIM(K).NE.TIMP .AND. K.GT.I1)THEN
     SUMT = 0.0
     SUMI = 0.0
    DO I=I1,K-1
     IF(SLA%FLG(I).EQ.1)THEN
      SUMT = SUMT + SLA%RES(I)
      SUMI = SUMI + 1.0
     ENDIF
    ENDDO
     IF(SUMI.GT.0.) SUMT = SUMT/SUMI
    DO I=I1,K-1
     SLA%RES(I) = SLA%RES(I) - SUMT
     SLA%BIA(I) = SUMT
    ENDDO
     TIMP = SLA%TIM(K)
     I1 = K
  ELSE IF(K.EQ.SLA%NO .AND. K.GE.I1)THEN
     SUMT = 0.0
     SUMI = 0.0
    DO I=I1,K
     IF(SLA%FLG(I).EQ.1)THEN
      SUMT = SUMT + SLA%RES(I)
      SUMI = SUMI + 1.0
     ENDIF
    ENDDO
     IF(SUMI.GT.0.) SUMT = SUMT/SUMI
    DO I=I1,K
     SLA%RES(I) = SLA%RES(I) - SUMT
     SLA%BIA(I) = SUMT
    ENDDO
  ENDIF
 ENDDO

! RESIDUAL CHECK
  DO K=1,SLA%NO
   IF(ABS(SLA%RES(K)).GT.0.3) SLA%FLG(K) = 0
  ENDDO

 ENDDO ! ITER

! ---
! COUNT GOOD OBSERVATIONS
    SLA%NC = 0
  DO K=1,SLA%NO
   IF(SLA%FLG(K).EQ.1)THEN
    SLA%NC = SLA%NC + 1
   ELSE
    SLA%FLC(K) = 0
    SLA%BIA(K) = 0.
    SLA%RES(K) = 0.
    SLA%INC(K) = 0.
    SLA%B_A(K) = 0.
    SLA%PQ1(K) = 0.
    SLA%PQ2(K) = 0.
    SLA%PQ3(K) = 0.
    SLA%PQ4(K) = 0.
   ENDIF
  ENDDO

  SLA%FLC(:) = SLA%FLG(:)

#endif
END SUBROUTINE GET_OBS_SLA
