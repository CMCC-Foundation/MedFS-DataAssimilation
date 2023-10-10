SUBROUTINE OBS_VEL

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR SLA                                 !
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

    I=VEL%IB(K)
    J=VEL%JB(K)

     VEL%INC(K) = VEL%PQ1(K) * GRD%UVL(I  ,J  ,VEL%KDP) +       &
                  VEL%PQ2(K) * GRD%UVL(I+1,J  ,VEL%KDP) +       &
                  VEL%PQ3(K) * GRD%UVL(I  ,J+1,VEL%KDP) +       &
                  VEL%PQ4(K) * GRD%UVL(I+1,J+1,VEL%KDP)


  ELSE IF(VEL%FLC(K).EQ.1 .AND. VEL%PAR(K).EQ.2)THEN

    I=VEL%IB(K)
    J=VEL%JB(K)

     VEL%INC(K) = VEL%PQ1(K) * GRD%VVL(I  ,J  ,VEL%KDP) +       &
                  VEL%PQ2(K) * GRD%VVL(I+1,J  ,VEL%KDP) +       &
                  VEL%PQ3(K) * GRD%VVL(I  ,J+1,VEL%KDP) +       &
                  VEL%PQ4(K) * GRD%VVL(I+1,J+1,VEL%KDP)


  ENDIF

 ENDDO

END SUBROUTINE OBS_VEL
