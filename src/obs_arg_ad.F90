SUBROUTINE OBS_ARG_AD

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR ARGO FLOATS - ADJOINT
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

 INTEGER(I4)   ::  I, J, K, KK

 DO KK = 1,ARG%NO

  IF(ARG%FLC(KK).EQ.1 .AND. ARG%PAR(KK).EQ.1 )THEN

    OBS%K = OBS%K + 1

    I=ARG%IB(KK)
    J=ARG%JB(KK)
    K=ARG%KB(KK)

    GRD%TEM_AD(I  ,J  ,K  ) = GRD%TEM_AD(I  ,J  ,K  ) + ARG%PQ1(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J  ,K  ) = GRD%TEM_AD(I+1,J  ,K  ) + ARG%PQ2(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I  ,J+1,K  ) = GRD%TEM_AD(I  ,J+1,K  ) + ARG%PQ3(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J+1,K  ) = GRD%TEM_AD(I+1,J+1,K  ) + ARG%PQ4(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I  ,J  ,K+1) = GRD%TEM_AD(I  ,J  ,K+1) + ARG%PQ5(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J  ,K+1) = GRD%TEM_AD(I+1,J  ,K+1) + ARG%PQ6(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I  ,J+1,K+1) = GRD%TEM_AD(I  ,J+1,K+1) + ARG%PQ7(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J+1,K+1) = GRD%TEM_AD(I+1,J+1,K+1) + ARG%PQ8(KK) * OBS%GRA(OBS%K)

  ELSE IF(ARG%FLC(KK).EQ.1 .AND. ARG%PAR(KK).EQ.2 )THEN

    OBS%K = OBS%K + 1

    I=ARG%IB(KK)
    J=ARG%JB(KK)
    K=ARG%KB(KK)

    GRD%SAL_AD(I  ,J  ,K  ) = GRD%SAL_AD(I  ,J  ,K  ) + ARG%PQ1(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I+1,J  ,K  ) = GRD%SAL_AD(I+1,J  ,K  ) + ARG%PQ2(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I  ,J+1,K  ) = GRD%SAL_AD(I  ,J+1,K  ) + ARG%PQ3(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I+1,J+1,K  ) = GRD%SAL_AD(I+1,J+1,K  ) + ARG%PQ4(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I  ,J  ,K+1) = GRD%SAL_AD(I  ,J  ,K+1) + ARG%PQ5(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I+1,J  ,K+1) = GRD%SAL_AD(I+1,J  ,K+1) + ARG%PQ6(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I  ,J+1,K+1) = GRD%SAL_AD(I  ,J+1,K+1) + ARG%PQ7(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I+1,J+1,K+1) = GRD%SAL_AD(I+1,J+1,K+1) + ARG%PQ8(KK) * OBS%GRA(OBS%K)

  ENDIF

 ENDDO


END SUBROUTINE OBS_ARG_AD
