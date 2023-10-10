SUBROUTINE OBS_GLD_AD

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR GLIDERS - ADJOINT
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

    OBS%K = OBS%K + 1

    I=GLD%IB(KK)
    J=GLD%JB(KK)
    K=GLD%KB(KK)

    GRD%TEM_AD(I  ,J  ,K  ) = GRD%TEM_AD(I  ,J  ,K  ) + GLD%PQ1(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J  ,K  ) = GRD%TEM_AD(I+1,J  ,K  ) + GLD%PQ2(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I  ,J+1,K  ) = GRD%TEM_AD(I  ,J+1,K  ) + GLD%PQ3(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J+1,K  ) = GRD%TEM_AD(I+1,J+1,K  ) + GLD%PQ4(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I  ,J  ,K+1) = GRD%TEM_AD(I  ,J  ,K+1) + GLD%PQ5(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J  ,K+1) = GRD%TEM_AD(I+1,J  ,K+1) + GLD%PQ6(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I  ,J+1,K+1) = GRD%TEM_AD(I  ,J+1,K+1) + GLD%PQ7(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J+1,K+1) = GRD%TEM_AD(I+1,J+1,K+1) + GLD%PQ8(KK) * OBS%GRA(OBS%K)

  ELSE IF(GLD%FLC(KK).EQ.1 .AND. GLD%PAR(KK).EQ.2 )THEN

    OBS%K = OBS%K + 1

    I=GLD%IB(KK)
    J=GLD%JB(KK)
    K=GLD%KB(KK)

    GRD%SAL_AD(I  ,J  ,K  ) = GRD%SAL_AD(I  ,J  ,K  ) + GLD%PQ1(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I+1,J  ,K  ) = GRD%SAL_AD(I+1,J  ,K  ) + GLD%PQ2(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I  ,J+1,K  ) = GRD%SAL_AD(I  ,J+1,K  ) + GLD%PQ3(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I+1,J+1,K  ) = GRD%SAL_AD(I+1,J+1,K  ) + GLD%PQ4(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I  ,J  ,K+1) = GRD%SAL_AD(I  ,J  ,K+1) + GLD%PQ5(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I+1,J  ,K+1) = GRD%SAL_AD(I+1,J  ,K+1) + GLD%PQ6(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I  ,J+1,K+1) = GRD%SAL_AD(I  ,J+1,K+1) + GLD%PQ7(KK) * OBS%GRA(OBS%K)
    GRD%SAL_AD(I+1,J+1,K+1) = GRD%SAL_AD(I+1,J+1,K+1) + GLD%PQ8(KK) * OBS%GRA(OBS%K)

  ENDIF

 ENDDO


END SUBROUTINE OBS_GLD_AD
