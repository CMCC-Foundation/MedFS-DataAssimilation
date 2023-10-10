SUBROUTINE OBS_GLD

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR GLIDERS                             !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

 INTEGER(I4)   ::  I, J, K, KK

 DO KK = 1,GLD%NO

  IF(GLD%FLC(KK).EQ.1 .AND. GLD%PAR(KK).EQ.1 )THEN

    I=GLD%IB(KK)
    J=GLD%JB(KK)
    K=GLD%KB(KK)

    GLD%INC(KK) = GLD%PQ1(KK) * GRD%TEM(I  ,J  ,K  ) +       &
                  GLD%PQ2(KK) * GRD%TEM(I+1,J  ,K  ) +       &
                  GLD%PQ3(KK) * GRD%TEM(I  ,J+1,K  ) +       &
                  GLD%PQ4(KK) * GRD%TEM(I+1,J+1,K  ) +       &
                  GLD%PQ5(KK) * GRD%TEM(I  ,J  ,K+1) +       &
                  GLD%PQ6(KK) * GRD%TEM(I+1,J  ,K+1) +       &
                  GLD%PQ7(KK) * GRD%TEM(I  ,J+1,K+1) +       &
                  GLD%PQ8(KK) * GRD%TEM(I+1,J+1,K+1)

  ELSE IF(GLD%FLC(KK).EQ.1 .AND. GLD%PAR(KK).EQ.2 )THEN

    I=GLD%IB(KK)
    J=GLD%JB(KK)
    K=GLD%KB(KK)

    GLD%INC(KK) = GLD%PQ1(KK) * GRD%SAL(I  ,J  ,K  ) +       &
                  GLD%PQ2(KK) * GRD%SAL(I+1,J  ,K  ) +       &
                  GLD%PQ3(KK) * GRD%SAL(I  ,J+1,K  ) +       &
                  GLD%PQ4(KK) * GRD%SAL(I+1,J+1,K  ) +       &
                  GLD%PQ5(KK) * GRD%SAL(I  ,J  ,K+1) +       &
                  GLD%PQ6(KK) * GRD%SAL(I+1,J  ,K+1) +       &
                  GLD%PQ7(KK) * GRD%SAL(I  ,J+1,K+1) +       &
                  GLD%PQ8(KK) * GRD%SAL(I+1,J+1,K+1)

  ENDIF

 ENDDO

END SUBROUTINE OBS_GLD
