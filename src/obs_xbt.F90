SUBROUTINE OBS_XBT

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR XBT PROFILES                        !
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

    I=XBT%IB(KK)
    J=XBT%JB(KK)
    K=XBT%KB(KK)

    XBT%INC(KK) = XBT%PQ1(KK) * GRD%TEM(I  ,J  ,K  ) +       &
                  XBT%PQ2(KK) * GRD%TEM(I+1,J  ,K  ) +       &
                  XBT%PQ3(KK) * GRD%TEM(I  ,J+1,K  ) +       &
                  XBT%PQ4(KK) * GRD%TEM(I+1,J+1,K  ) +       &
                  XBT%PQ5(KK) * GRD%TEM(I  ,J  ,K+1) +       &
                  XBT%PQ6(KK) * GRD%TEM(I+1,J  ,K+1) +       &
                  XBT%PQ7(KK) * GRD%TEM(I  ,J+1,K+1) +       &
                  XBT%PQ8(KK) * GRD%TEM(I+1,J+1,K+1)

  ENDIF

 ENDDO

END SUBROUTINE OBS_XBT
