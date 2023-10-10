SUBROUTINE INVRT( IM, JM, FLD, MSK, RGH, A1, A2, A3, A4, A0,           &
                  BNM, OVR, RESEM, NCNT, ITR )

!-----------------------------------------------------------------------
!                                                                      !
! IMPLICIT SOLVER - OVERRELAXATION                                     !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND

  IMPLICIT NONE

 INTEGER(I4)    :: IM, JM

 REAL(R8)       :: FLD(IM,JM), MSK(IM,JM)
 REAL(R8)       :: A0(IM,JM), A1(IM,JM), A2(IM,JM), A3(IM,JM), A4(IM,JM)
 REAL(R8)       :: RGH(IM,JM)
 REAL(R8)       :: RES(IM,JM)

 REAL(R8)       :: BNM, RESER, OVR, RESEM
 INTEGER(I4)    :: I, J, ICNT, NCNT, ITR




            RESER = 1.E20

       ITR = 0

      DO ICNT=1,NCNT

       IF(RESER.GT.RESEM)THEN

         ITR = ITR + 1

         DO J=2,JM-1,2
          DO I=2,IM-1,2
            RES(I,J) = A1(I,J)*FLD(I+1,J  )+A2(I,J)*FLD(I-1,J  )+        &
                       A3(I,J)*FLD(I  ,J+1)+A4(I,J)*FLD(I  ,J-1)-        &
                       A0(I,J)*FLD(I  ,J  ) - RGH(I,J)
               FLD(I,J) = FLD(I,J)+OVR*RES(I,J)/A0(I,J)*MSK(I,J)
          ENDDO
         ENDDO
         DO J=3,JM-1,2
          DO I=3,IM-1,2
            RES(I,J) = A1(I,J)*FLD(I+1,J  )+A2(I,J)*FLD(I-1,J  )+        &
                       A3(I,J)*FLD(I  ,J+1)+A4(I,J)*FLD(I  ,J-1)-        &
                       A0(I,J)*FLD(I  ,J  ) - RGH(I,J)
               FLD(I,J) = FLD(I,J)+OVR*RES(I,J)/A0(I,J)*MSK(I,J)
          ENDDO
         ENDDO
         DO J=3,JM-1,2
          DO I=2,IM-1,2
            RES(I,J) = A1(I,J)*FLD(I+1,J  )+A2(I,J)*FLD(I-1,J  )+        &
                       A3(I,J)*FLD(I  ,J+1)+A4(I,J)*FLD(I  ,J-1)-        &
                       A0(I,J)*FLD(I  ,J  ) - RGH(I,J)
               FLD(I,J) = FLD(I,J)+OVR*RES(I,J)/A0(I,J)*MSK(I,J)
          ENDDO
         ENDDO
         DO J=2,JM-1,2
          DO I=3,IM-1,2
            RES(I,J) = A1(I,J)*FLD(I+1,J  )+A2(I,J)*FLD(I-1,J  )+        &
                       A3(I,J)*FLD(I  ,J+1)+A4(I,J)*FLD(I  ,J-1)-        &
                       A0(I,J)*FLD(I  ,J  ) - RGH(I,J)
               FLD(I,J) = FLD(I,J)+OVR*RES(I,J)/A0(I,J)*MSK(I,J)
          ENDDO
         ENDDO

            RESER = 0.0
         DO J=2,JM-1
          DO I=2,IM-1
            RESER = RESER + ABS(RES(I,J))
          ENDDO
         ENDDO
            RESER = RESER/BNM

       ENDIF

      ENDDO



END SUBROUTINE INVRT
