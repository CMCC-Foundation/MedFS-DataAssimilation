SUBROUTINE OBS_ARG

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR ARGO FLOATS                         !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR
 USE CTL_STR

 IMPLICIT NONE

 INTEGER(I4)   ::  I, J, K, KK

 DO KK = 1,ARG%NO

  IF(ARG%FLC(KK).EQ.1 .AND. ARG%PAR(KK).EQ.1 )THEN

    I=ARG%IB(KK)
    J=ARG%JB(KK)
    K=ARG%KB(KK)

    ARG%INC(KK) = ARG%PQ1(KK) * GRD%TEM(I  ,J  ,K  ) +       &
                  ARG%PQ2(KK) * GRD%TEM(I+1,J  ,K  ) +       &
                  ARG%PQ3(KK) * GRD%TEM(I  ,J+1,K  ) +       &
                  ARG%PQ4(KK) * GRD%TEM(I+1,J+1,K  ) +       &
                  ARG%PQ5(KK) * GRD%TEM(I  ,J  ,K+1) +       &
                  ARG%PQ6(KK) * GRD%TEM(I+1,J  ,K+1) +       &
                  ARG%PQ7(KK) * GRD%TEM(I  ,J+1,K+1) +       &
                  ARG%PQ8(KK) * GRD%TEM(I+1,J+1,K+1)

  ELSE IF(ARG%FLC(KK).EQ.1 .AND. ARG%PAR(KK).EQ.2 )THEN

    I=ARG%IB(KK)
    J=ARG%JB(KK)
    K=ARG%KB(KK)

    ARG%INC(KK) = ARG%PQ1(KK) * GRD%SAL(I  ,J  ,K  ) +       &
                  ARG%PQ2(KK) * GRD%SAL(I+1,J  ,K  ) +       &
                  ARG%PQ3(KK) * GRD%SAL(I  ,J+1,K  ) +       &
                  ARG%PQ4(KK) * GRD%SAL(I+1,J+1,K  ) +       &
                  ARG%PQ5(KK) * GRD%SAL(I  ,J  ,K+1) +       &
                  ARG%PQ6(KK) * GRD%SAL(I+1,J  ,K+1) +       &
                  ARG%PQ7(KK) * GRD%SAL(I  ,J+1,K+1) +       &
                  ARG%PQ8(KK) * GRD%SAL(I+1,J+1,K+1)

  ENDIF

 ENDDO

END SUBROUTINE OBS_ARG
