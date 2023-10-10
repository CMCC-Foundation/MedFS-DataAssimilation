SUBROUTINE CNV_INN

!-----------------------------------------------------------------------
!                                                                      !
! CONVERT W TO CORRECTION IN PHYSICAL SPACE                            !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

 USE SET_KND
 USE OBS_STR
 USE GRD_STR
 USE EOF_STR
 USE CTL_STR
 USE DRV_STR
 USE RUN
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL
 USE TLAD_VARS
 USE TLAD
 USE HYBRID
 USE CERES

 IMPLICIT NONE

 INTEGER(I4)    :: I,J,K, KK
 REAL(R8)    :: TMN,SMN, SPD,ONEN, ZA1, ZA2

 CALL MYFRTPROF_WALL('CNV_INN: CONVERT V INTO INNOVATIONS',0)

 DRV%DDA(DRV%KTR) = DRV%DDI(DRV%KTR)

 IF(LL_DIAHYBRID .AND. LL_HYBRID) LL_DIAH = .TRUE.

! --------
! CONVERT THE CONTROL VECTOR TO V
   CALL CNV_CTV

   IF( NCONF .EQ. 202 ) THEN
     CALL VER_HOR202
   ELSE
     IF(LL_HYBRID_CR) CALL VER_HORF
     CALL VER_HOR
   ENDIF

   IF( LL_HYBRID .AND.  LL_ALPHA_CTL ) THEN
      IF(LL_PRECON_ALPHA) THEN
        ZA1 = ALPHAB(1) + ALPHA(1)*SQRT(ZALPHA_VAR)
      ELSE
        ZA1 = ALPHAB(1) + ALPHA(1)
      ENDIF
      ZA2 = SQRT(1._R8-ZA1*ZA1)
      OPEN(795,FILE='ALPHA_CTL.dat',STATUS='NEW')
      WRITE(795,'(4F10.6)') ALPHAB(1),ZA1,ZA1**2,ZA2**2
      CLOSE(795)
   ENDIF


! --------
! SIMPLIFIED TL MODEL
  IF ( LL_WEAKLY ) CALL CTRL2ERR_TL(1)

  IF ( LL_4DVAR .OR. LL_TLMODEL ) CALL TLMOD(1)

  IF(LL_CERES) CALL CERES_OO_TLF(GRD%IM,GRD%JM,GRD%KM,GRD%TEM)

  IF ( LL_TLMODEL .AND. LL_ADMODEL ) THEN
     GRD%TEM_AD=GRD%TEM
     GRD%SAL_AD=GRD%SAL
     CALL ADMOD(1)
  ENDIF

 CALL MYFRTPROF_WALL('CNV_INN: CONVERT V INTO INNOVATIONS',1)
END SUBROUTINE CNV_INN
