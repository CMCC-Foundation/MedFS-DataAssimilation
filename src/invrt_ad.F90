SUBROUTINE INVRT_AD( IM, JM, FLD, MSK, RGH, A1, A2, A3, A4, A0,           &
                  BNM, OVR, RESEM, NCNT, ITR )

!-----------------------------------------------------------------------
!                                                                      !
! IMPLICIT SOLVER - OVERRELAXATION (ADJOINT)                           !
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


   RES(:,:) = 0.0

      DO ICNT=ITR,1,-1

         DO J=2,JM-1,2
          DO I=3,IM-1,2
            RES(I,J) = RES(I,J) + FLD(I,J)*OVR/A0(I,J)*MSK(I,J)
            FLD(I,J-1) = FLD(I,J-1) + RES(I,J)*A4(I,J)
            FLD(I,J+1) = FLD(I,J+1) + RES(I,J)*A3(I,J)
            FLD(I-1,J) = FLD(I-1,J) + RES(I,J)*A2(I,J)
            FLD(I+1,J) = FLD(I+1,J) + RES(I,J)*A1(I,J)
            FLD(I,J) = FLD(I,J) - RES(I,J)*A0(I,J)
            RGH(I,J) = RGH(I,J) - RES(I,J)
            RES(I,J) = 0.0
          ENDDO
         ENDDO
         DO J=3,JM-1,2
          DO I=2,IM-1,2
            RES(I,J) = RES(I,J) + FLD(I,J)*OVR/A0(I,J)*MSK(I,J)
            FLD(I,J-1) = FLD(I,J-1) + RES(I,J)*A4(I,J)
            FLD(I,J+1) = FLD(I,J+1) + RES(I,J)*A3(I,J)
            FLD(I-1,J) = FLD(I-1,J) + RES(I,J)*A2(I,J)
            FLD(I+1,J) = FLD(I+1,J) + RES(I,J)*A1(I,J)
            FLD(I,J) = FLD(I,J) - RES(I,J)*A0(I,J)
            RGH(I,J) = RGH(I,J) - RES(I,J)
            RES(I,J) = 0.0
          ENDDO
         ENDDO
         DO J=3,JM-1,2
          DO I=3,IM-1,2
            RES(I,J) = RES(I,J) + FLD(I,J)*OVR/A0(I,J)*MSK(I,J)
            FLD(I,J-1) = FLD(I,J-1) + RES(I,J)*A4(I,J)
            FLD(I,J+1) = FLD(I,J+1) + RES(I,J)*A3(I,J)
            FLD(I-1,J) = FLD(I-1,J) + RES(I,J)*A2(I,J)
            FLD(I+1,J) = FLD(I+1,J) + RES(I,J)*A1(I,J)
            FLD(I,J) = FLD(I,J) - RES(I,J)*A0(I,J)
            RGH(I,J) = RGH(I,J) - RES(I,J)
            RES(I,J) = 0.0
          ENDDO
         ENDDO
         DO J=2,JM-1,2
          DO I=2,IM-1,2
            RES(I,J) = RES(I,J) + FLD(I,J)*OVR/A0(I,J)*MSK(I,J)
            FLD(I,J-1) = FLD(I,J-1) + RES(I,J)*A4(I,J)
            FLD(I,J+1) = FLD(I,J+1) + RES(I,J)*A3(I,J)
            FLD(I-1,J) = FLD(I-1,J) + RES(I,J)*A2(I,J)
            FLD(I+1,J) = FLD(I+1,J) + RES(I,J)*A1(I,J)
            FLD(I,J) = FLD(I,J) - RES(I,J)*A0(I,J)
            RGH(I,J) = RGH(I,J) - RES(I,J)
            RES(I,J) = 0.0
          ENDDO
         ENDDO

            FLD(:,:) = FLD(:,:) * MSK(:,:)

      ENDDO



END SUBROUTINE INVRT_AD
