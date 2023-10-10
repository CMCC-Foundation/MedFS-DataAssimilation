SUBROUTINE OBS_TRJ_AD

!-----------------------------------------------------------------------
!                                                                      !
! CALL TRAJECTORY MODEL                                                !
!                                                                      !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

 INTEGER(I4)   ::  K

 DO K=1,TRJ%JPN

  IF(TRJ%FLC(K).EQ.1 .AND. TRJ%PAR(K).EQ.1)THEN

     OBS%K = OBS%K + 1

     TRJ%XTL_AD(K) = OBS%GRA(OBS%K)

  ENDIF

 ENDDO

 DO K=1,TRJ%JPN

  IF(TRJ%FLC(TRJ%JPN+K).EQ.1 .AND. TRJ%PAR(TRJ%JPN+K).EQ.2)THEN

     OBS%K = OBS%K + 1

     TRJ%YTL_AD(K) = OBS%GRA(OBS%K)

  ENDIF

 ENDDO


  PRINT*,'INPUT FOR ADJOINT: ', TRJ%XTL_AD, TRJ%YTL_AD

  CALL MOD_TRJ_AD( GRD%IM,GRD%JM,TRJ%UMN,TRJ%VMN,GRD%DX,GRD%DY,      &
                   TRJ%JPT,TRJ%JPN,TRJ%XMN,TRJ%YMN,TRJ%TIM,          &
                   GRD%UVL_AD(1,1,TRJ%KDP),GRD%VVL_AD(1,1,TRJ%KDP),  &
                   TRJ%XTL_AD,TRJ%YTL_AD )
  PRINT*,'INPUT FOR ADJOINT: ', TRJ%XTL_AD, TRJ%YTL_AD


!   OPEN(101,FILE='CORRECTIONS_1.DAT',FORM='UNFORMATTED')
!    WRITE(101) GRD%ETA
!   DO K=1,GRD%KM
!    WRITE(101) GRD%TEM(:,:,K)
!   ENDDO
!   DO K=1,GRD%KM
!    WRITE(101) GRD%SAL(:,:,K)
!   ENDDO
!   DO K=1,GRD%KM
!    WRITE(101) GRD%UVL_AD(:,:,TRJ%KDP)
!   ENDDO
!   DO K=1,GRD%KM
!    WRITE(101) GRD%VVL_AD(:,:,TRJ%KDP)
!   ENDDO
!   DO K=1,GRD%KM
!!    WRITE(101) REAL(GRD%MSK(:,:,K))
!   ENDDO
!   CLOSE(101)
!   STOP

END SUBROUTINE OBS_TRJ_AD
