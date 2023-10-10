SUBROUTINE OBS_XBT_AD

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

 DO KK = 1,XBT%NO

  IF(XBT%FLC(KK).EQ.1 .AND. XBT%PAR(KK).EQ.1 )THEN

    OBS%K = OBS%K + 1

    I=XBT%IB(KK)
    J=XBT%JB(KK)
    K=XBT%KB(KK)

    GRD%TEM_AD(I  ,J  ,K  ) = GRD%TEM_AD(I  ,J  ,K  ) + XBT%PQ1(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J  ,K  ) = GRD%TEM_AD(I+1,J  ,K  ) + XBT%PQ2(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I  ,J+1,K  ) = GRD%TEM_AD(I  ,J+1,K  ) + XBT%PQ3(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J+1,K  ) = GRD%TEM_AD(I+1,J+1,K  ) + XBT%PQ4(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I  ,J  ,K+1) = GRD%TEM_AD(I  ,J  ,K+1) + XBT%PQ5(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J  ,K+1) = GRD%TEM_AD(I+1,J  ,K+1) + XBT%PQ6(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I  ,J+1,K+1) = GRD%TEM_AD(I  ,J+1,K+1) + XBT%PQ7(KK) * OBS%GRA(OBS%K)
    GRD%TEM_AD(I+1,J+1,K+1) = GRD%TEM_AD(I+1,J+1,K+1) + XBT%PQ8(KK) * OBS%GRA(OBS%K)


  ENDIF

 ENDDO


END SUBROUTINE OBS_XBT_AD
