SUBROUTINE OBS_TRJ

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

  PRINT*,'INPUT:'
  PRINT*,GRD%IM,GRD%JM,TRJ%JPT,TRJ%JPN,TRJ%XMN,TRJ%YMN,TRJ%TIM
  CALL MOD_TRJ_TL( GRD%IM,GRD%JM,TRJ%UMN,TRJ%VMN,GRD%DX,GRD%DY,      &
                   TRJ%JPT,TRJ%JPN,TRJ%XMN,TRJ%YMN,TRJ%TIM,          &
                   GRD%UVL(1,1,TRJ%KDP),GRD%VVL(1,1,TRJ%KDP),TRJ%XTL,TRJ%YTL )
  PRINT*,'INCREMENTS: ',TRJ%XTL,TRJ%YTL


 DO K=1,TRJ%JPN

  IF(TRJ%FLC(K).EQ.1 .AND. TRJ%PAR(K).EQ.1)THEN

     TRJ%INC(K) = TRJ%XTL(K)

  ENDIF

 ENDDO

 DO K=1,TRJ%JPN

  IF(TRJ%FLC(TRJ%JPN+K).EQ.1 .AND. TRJ%PAR(TRJ%JPN+K).EQ.2)THEN

     TRJ%INC(TRJ%JPN+K) = TRJ%YTL(K)

  ENDIF

 ENDDO

END SUBROUTINE OBS_TRJ
