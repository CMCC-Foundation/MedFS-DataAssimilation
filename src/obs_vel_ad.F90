SUBROUTINE OBS_VEL_AD

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR SLA - ADJOINT                       !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

 INTEGER(I4)   ::  I, J, K

 DO K=1,VEL%NO

  IF(VEL%FLC(K).EQ.1 .AND. VEL%PAR(K).EQ.1)THEN

    OBS%K = OBS%K + 1

    I=VEL%IB(K)
    J=VEL%JB(K)

    GRD%UVL_AD(I  ,J  ,VEL%KDP) = GRD%UVL_AD(I  ,J  ,VEL%KDP) + VEL%PQ1(K) * OBS%GRA(OBS%K)
    GRD%UVL_AD(I+1,J  ,VEL%KDP) = GRD%UVL_AD(I+1,J  ,VEL%KDP) + VEL%PQ2(K) * OBS%GRA(OBS%K)
    GRD%UVL_AD(I  ,J+1,VEL%KDP) = GRD%UVL_AD(I  ,J+1,VEL%KDP) + VEL%PQ3(K) * OBS%GRA(OBS%K)
    GRD%UVL_AD(I+1,J+1,VEL%KDP) = GRD%UVL_AD(I+1,J+1,VEL%KDP) + VEL%PQ4(K) * OBS%GRA(OBS%K)

  ELSE IF(VEL%FLC(K).EQ.1 .AND. VEL%PAR(K).EQ.2)THEN

    OBS%K = OBS%K + 1

    I=VEL%IB(K)
    J=VEL%JB(K)

    GRD%VVL_AD(I  ,J  ,VEL%KDP) = GRD%VVL_AD(I  ,J  ,VEL%KDP) + VEL%PQ1(K) * OBS%GRA(OBS%K)
    GRD%VVL_AD(I+1,J  ,VEL%KDP) = GRD%VVL_AD(I+1,J  ,VEL%KDP) + VEL%PQ2(K) * OBS%GRA(OBS%K)
    GRD%VVL_AD(I  ,J+1,VEL%KDP) = GRD%VVL_AD(I  ,J+1,VEL%KDP) + VEL%PQ3(K) * OBS%GRA(OBS%K)
    GRD%VVL_AD(I+1,J+1,VEL%KDP) = GRD%VVL_AD(I+1,J+1,VEL%KDP) + VEL%PQ4(K) * OBS%GRA(OBS%K)

  ENDIF

 ENDDO

END SUBROUTINE OBS_VEL_AD
