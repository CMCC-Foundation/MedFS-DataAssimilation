SUBROUTINE BAR_MOD_AD

!-----------------------------------------------------------------------
!                                                                      !
! BAROTROPIC MODEL  (ADJOINT)                                          !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR

  IMPLICIT NONE

  INTEGER(I4)    :: I, J, KSTP


                 BMD%BXBY(:,:) = 0.
                 BMD%CU(:,:) = 0.
                 BMD%CV(:,:) = 0.
                 BMD%DIV(:,:) = 0.
                 BMD%DUX(:,:) = 0.
                 BMD%DUY(:,:) = 0.
                 BMD%DVX(:,:) = 0.
                 BMD%DVY(:,:) = 0.
                 BMD%ETA(:,:) = 0.
                 BMD%ETB(:,:) = 0.
                 BMD%ETX(:,:) = 0.
                 BMD%ETY(:,:) = 0.
                 BMD%RGH(:,:) = 0.
                 BMD%UA(:,:) = 0.
                 BMD%UB(:,:) = 0.
                 BMD%UN(:,:) = 0.
                 BMD%VA(:,:) = 0.
                 BMD%VB(:,:) = 0.
                 BMD%VN(:,:) = 0.
                 GRD%BX(:,:) = 0.
                 GRD%BY(:,:) = 0.

!---

     DO KSTP = BMD%NSTPS, 1, -1



        IF (MOD(KSTP,BMD%NSTPA) .EQ. 0) THEN

          BMD%ETM(:,:) = 0.0
           BMD%VM(:,:) = 0.0
           BMD%UM(:,:) = 0.0

           BMD%ETM(:,:) = GRD%ETA_AD(:,:)
           GRD%ETA_AD(:,:) = 0.0
 
           BMD%VM(:,:) =  BMD%VM(:,:)/BMD%NSTPA
           BMD%UM(:,:) =  BMD%UM(:,:)/BMD%NSTPA
          BMD%ETM(:,:) = BMD%ETM(:,:)/BMD%NSTPA

        ENDIF

          BMD%ETB(:,:) = BMD%ETB(:,:) + BMD%ETM(:,:)
           BMD%UB(:,:) =  BMD%UB(:,:) +  BMD%UM(:,:)
           BMD%VB(:,:) =  BMD%VB(:,:) +  BMD%VM(:,:)

           BMD%UA(:,:) =  BMD%UA(:,:) +  BMD%UN(:,:)
           BMD%VA(:,:) =  BMD%VA(:,:) +  BMD%VN(:,:)
           BMD%UN(:,:) = 0.0
           BMD%VN(:,:) = 0.0

          BMD%ETA(:,:) = BMD%ETA(:,:) + BMD%ETB(:,:)
           BMD%UN(:,:) =  BMD%UN(:,:) +  BMD%UB(:,:)
           BMD%VN(:,:) =  BMD%VN(:,:) +  BMD%VB(:,:)
          BMD%ETB(:,:) = 0.0
           BMD%UB(:,:) = 0.0
           BMD%VB(:,:) = 0.0

           BMD%UA(:,:) =  BMD%UA(:,:) + BMD%UN(:,:)*0.05
           BMD%UB(:,:) =  BMD%UB(:,:) + BMD%UN(:,:)*0.05
           BMD%UN(:,:) =  BMD%UN(:,:) * (1.0 - 2.0*0.05)
           BMD%VA(:,:) =  BMD%VA(:,:) + BMD%VN(:,:)*0.05
           BMD%VB(:,:) =  BMD%VB(:,:) + BMD%VN(:,:)*0.05
           BMD%VN(:,:) =  BMD%VN(:,:) * (1.0 - 2.0*0.05)

       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          BMD%DUX(I+1,J) = BMD%DUX(I+1,J) + BMD%DF2*BMD%UA(I,J)/BMD%DXU(I,J)*BMD%MSU(I,J)
          BMD%DUX(I  ,J) = BMD%DUX(I  ,J) - BMD%DF2*BMD%UA(I,J)/BMD%DXU(I,J)*BMD%MSU(I,J)
          BMD%DUY(I,J+1) = BMD%DUY(I,J+1) + BMD%DF2*BMD%UA(I,J)/BMD%DYU(I,J)*BMD%MSU(I,J)
          BMD%DUY(I,J  ) = BMD%DUY(I,J  ) - BMD%DF2*BMD%UA(I,J)/BMD%DYU(I,J)*BMD%MSU(I,J)
        ENDDO
       ENDDO
       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          BMD%DVX(I+1,J) = BMD%DVX(I+1,J) + BMD%DF2*BMD%VA(I,J)/BMD%DXV(I,J)*BMD%MSV(I,J)
          BMD%DVX(I  ,J) = BMD%DVX(I  ,J) - BMD%DF2*BMD%VA(I,J)/BMD%DXV(I,J)*BMD%MSV(I,J)
          BMD%DVY(I,J+1) = BMD%DVY(I,J+1) + BMD%DF2*BMD%VA(I,J)/BMD%DYV(I,J)*BMD%MSV(I,J)
          BMD%DVY(I,J  ) = BMD%DVY(I,J  ) - BMD%DF2*BMD%VA(I,J)/BMD%DYV(I,J)*BMD%MSV(I,J)
        ENDDO
       ENDDO

        DO J=1,GRD%JM
         DO I=2,GRD%IM
          BMD%UA(I  ,J) = BMD%UA(I  ,J) + BMD%DUX(I,J)/BMD%DXU(I,J)
          BMD%UA(I-1,J) = BMD%UA(I-1,J) - BMD%DUX(I,J)/BMD%DXU(I,J)
          BMD%VA(I  ,J) = BMD%VA(I  ,J) + BMD%DVX(I,J)/BMD%DXV(I,J)
          BMD%VA(I-1,J) = BMD%VA(I-1,J) - BMD%DVX(I,J)/BMD%DXV(I,J)
          BMD%DUX(I,J) = 0.0
          BMD%DVX(I,J) = 0.0
         ENDDO
        ENDDO
        DO J=2,GRD%JM
         DO I=1,GRD%IM
          BMD%UA(I,J  ) = BMD%UA(I,J  ) + BMD%DUY(I,J)/BMD%DYU(I,J)
          BMD%UA(I,J-1) = BMD%UA(I,J-1) - BMD%DUY(I,J)/BMD%DYU(I,J)
          BMD%VA(I,J  ) = BMD%VA(I,J  ) + BMD%DVY(I,J)/BMD%DYV(I,J)
          BMD%VA(I,J-1) = BMD%VA(I,J-1) - BMD%DVY(I,J)/BMD%DYV(I,J)
          BMD%DUY(I,J) = 0.0
          BMD%DVY(I,J) = 0.0
         ENDDO
        ENDDO


       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          GRD%BX(I,J)    = GRD%BX(I,J) - BMD%DT*BMD%UA(I,J)*BMD%MSU(I,J)
          BMD%CU(I,J)    = BMD%CU(I,J) - BMD%DT*BMD%UA(I,J)*BMD%MSU(I,J)
          BMD%DUX(I+1,J) = BMD%DUX(I+1,J) + BMD%DF1*BMD%UA(I,J)/BMD%DXU(I,J)*BMD%MSU(I,J)
          BMD%DUX(I  ,J) = BMD%DUX(I  ,J) - BMD%DF1*BMD%UA(I,J)/BMD%DXU(I,J)*BMD%MSU(I,J)
          BMD%DUY(I,J+1) = BMD%DUY(I,J+1) + BMD%DF1*BMD%UA(I,J)/BMD%DYU(I,J)*BMD%MSU(I,J)
          BMD%DUY(I,J  ) = BMD%DUY(I,J  ) - BMD%DF1*BMD%UA(I,J)/BMD%DYU(I,J)*BMD%MSU(I,J)
          BMD%ETA(I  ,J) = BMD%ETA(I  ,J) - BMD%DT*BMD%ALP1*BMD%G*BMD%HGU(I,J)*BMD%UA(I,J)/BMD%DXU(I,J)*BMD%MSU(I,J)
          BMD%ETA(I-1,J) = BMD%ETA(I-1,J) + BMD%DT*BMD%ALP1*BMD%G*BMD%HGU(I,J)*BMD%UA(I,J)/BMD%DXU(I,J)*BMD%MSU(I,J)
          BMD%ETX(I,J)   = BMD%ETX(I,J) - BMD%UA(I,J)
          BMD%UB(I,J)    = BMD%UB(I,J) + BMD%UA(I,J)
          BMD%UA(I,J)    = 0.0
        ENDDO
       ENDDO

       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          GRD%BY(I,J)    = GRD%BY(I,J) - BMD%DT*BMD%VA(I,J)*BMD%MSV(I,J)
          BMD%CV(I,J)    = BMD%CV(I,J) - BMD%DT*BMD%VA(I,J)*BMD%MSV(I,J)
          BMD%DVX(I+1,J) = BMD%DVX(I+1,J) + BMD%DF1*BMD%VA(I,J)/BMD%DXV(I,J)*BMD%MSV(I,J)
          BMD%DVX(I  ,J) = BMD%DVX(I  ,J) - BMD%DF1*BMD%VA(I,J)/BMD%DXV(I,J)*BMD%MSV(I,J)
          BMD%DVY(I,J+1) = BMD%DVY(I,J+1) + BMD%DF1*BMD%VA(I,J)/BMD%DYV(I,J)*BMD%MSV(I,J)
          BMD%DVY(I,J  ) = BMD%DVY(I,J  ) - BMD%DF1*BMD%VA(I,J)/BMD%DYV(I,J)*BMD%MSV(I,J)
          BMD%ETA(I,J  ) = BMD%ETA(I,J  ) - BMD%DT*BMD%ALP1*BMD%G*BMD%HGV(I,J)*BMD%VA(I,J)/BMD%DYV(I,J)*BMD%MSV(I,J)
          BMD%ETA(I,J-1) = BMD%ETA(I,J-1) + BMD%DT*BMD%ALP1*BMD%G*BMD%HGV(I,J)*BMD%VA(I,J)/BMD%DYV(I,J)*BMD%MSV(I,J)
          BMD%ETY(I,J)   = BMD%ETY(I,J) - BMD%VA(I,J)
          BMD%VB(I,J)    = BMD%VB(I,J) + BMD%VA(I,J)
          BMD%VA(I,J)    = 0.0
        ENDDO
       ENDDO

          BMD%UB(2:GRD%IM  ,:) = BMD%UB(2:GRD%IM  ,:) + BMD%DUX(2:GRD%IM,:)/BMD%DXU(2:GRD%IM,:)
          BMD%UB(1:GRD%IM-1,:) = BMD%UB(1:GRD%IM-1,:) - BMD%DUX(2:GRD%IM,:)/BMD%DXU(2:GRD%IM,:)
          BMD%DUX(2:GRD%IM,:)  = 0.0
          BMD%VB(2:GRD%IM  ,:) = BMD%VB(2:GRD%IM  ,:) + BMD%DVX(2:GRD%IM,:)/BMD%DXV(2:GRD%IM,:)
          BMD%VB(1:GRD%IM-1,:) = BMD%VB(1:GRD%IM-1,:) - BMD%DVX(2:GRD%IM,:)/BMD%DXV(2:GRD%IM,:)
          BMD%DVX(2:GRD%IM,:)  = 0.0

          BMD%UB(:,2:GRD%JM  ) = BMD%UB(:,2:GRD%JM  ) + BMD%DUY(:,2:GRD%JM)/BMD%DYU(:,2:GRD%JM)
          BMD%UB(:,1:GRD%JM-1) = BMD%UB(:,1:GRD%JM-1) - BMD%DUY(:,2:GRD%JM)/BMD%DYU(:,2:GRD%JM)
          BMD%DUY(2:GRD%IM,:)  = 0.0
          BMD%VB(:,2:GRD%JM  ) = BMD%VB(:,2:GRD%JM  ) + BMD%DVY(:,2:GRD%JM)/BMD%DYV(:,2:GRD%JM)
          BMD%VB(:,1:GRD%JM-1) = BMD%VB(:,1:GRD%JM-1) - BMD%DVY(:,2:GRD%JM)/BMD%DYV(:,2:GRD%JM)
          BMD%DVY(2:GRD%IM,:)  = 0.0

       DO J=1,GRD%JM-1
        DO I=2,GRD%IM-1
         BMD%VN(I-1,J+1) = BMD%VN(I-1,J+1) - BMD%CU(I,J)*BMD%DXV(I-1,J+1)*GRD%F(I,J+1)/BMD%DXU(I,J)*0.25
         BMD%VN(I-1,J  ) = BMD%VN(I-1,J  ) - BMD%CU(I,J)*BMD%DXV(I-1,J  )*GRD%F(I,J  )/BMD%DXU(I,J)*0.25
         BMD%VN(I  ,J+1) = BMD%VN(I  ,J+1) - BMD%CU(I,J)*BMD%DXV(I  ,J+1)*GRD%F(I,J+1)/BMD%DXU(I,J)*0.25
         BMD%VN(I  ,J  ) = BMD%VN(I  ,J  ) - BMD%CU(I,J)*BMD%DXV(I  ,J  )*GRD%F(I,J  )/BMD%DXU(I,J)*0.25
         BMD%CU(I,J)     = 0.0
        ENDDO
       ENDDO
       DO J=2,GRD%JM-1
        DO I=1,GRD%IM-1
         BMD%UN(I+1,J-1) = BMD%UN(I+1,J-1) + BMD%CV(I,J)*BMD%DYU(I+1,J-1)*GRD%F(I+1,J)/BMD%DYV(I,J)*0.25
         BMD%UN(I  ,J-1) = BMD%UN(I  ,J-1) + BMD%CV(I,J)*BMD%DYU(I  ,J-1)*GRD%F(I  ,J)/BMD%DYV(I,J)*0.25
         BMD%UN(I+1,J  ) = BMD%UN(I+1,J  ) + BMD%CV(I,J)*BMD%DYU(I+1,J  )*GRD%F(I+1,J)/BMD%DYV(I,J)*0.25
         BMD%UN(I  ,J  ) = BMD%UN(I  ,J  ) + BMD%CV(I,J)*BMD%DYU(I  ,J  )*GRD%F(I  ,J)/BMD%DYV(I,J)*0.25
         BMD%CV(I,J)     = 0.0
        ENDDO
       ENDDO

       CALL INVRT_AD( GRD%IM, GRD%JM, BMD%ETA, BMD%MST, BMD%RGH, BMD%A1, BMD%A2, BMD%A3, BMD%A4, BMD%A0,           &
                      BMD%BNM, BMD%OVR, BMD%RESEM, BMD%NCNT, BMD%ITR(KSTP) )

       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
         BMD%BXBY(I,J)  = BMD%BXBY(I,J) + BMD%RGH(I,J)
         BMD%DIV(I,J)   = BMD%DIV(I,J) + BMD%DT*BMD%RGH(I,J)
         BMD%ETB(I,J)   = BMD%ETB(I,J) - BMD%RGH(I,J)
         BMD%ETX(I+1,J) = BMD%ETX(I+1,J) - BMD%ALP1*BMD%DT*BMD%RGH(I,J)/GRD%DX(I,J)
         BMD%ETX(I  ,J) = BMD%ETX(I  ,J) + BMD%ALP1*BMD%DT*BMD%RGH(I,J)/GRD%DX(I,J)
         BMD%ETY(I,J+1) = BMD%ETY(I,J+1) - BMD%ALP1*BMD%DT*BMD%RGH(I,J)/GRD%DY(I,J)
         BMD%ETY(I,J  ) = BMD%ETY(I,J  ) + BMD%ALP1*BMD%DT*BMD%RGH(I,J)/GRD%DY(I,J)
         BMD%RGH(I,J  ) = 0.0
        ENDDO
       ENDDO

       DO J= 2,GRD%JM
        DO I= 2,GRD%IM-1
          BMD%ETB(I,J  ) = BMD%ETB(I,J  ) + BMD%ALP2*BMD%DT*BMD%G*BMD%HGV(I,J)*BMD%ETY(I,J)/BMD%DYV(I,J) * BMD%MSV(I,J)
          BMD%ETB(I,J-1) = BMD%ETB(I,J-1) - BMD%ALP2*BMD%DT*BMD%G*BMD%HGV(I,J)*BMD%ETY(I,J)/BMD%DYV(I,J) * BMD%MSV(I,J)
          BMD%ETY(I,J)   = 0.0
        ENDDO
       ENDDO
       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM
          BMD%ETB(I  ,J) = BMD%ETB(I  ,J) + BMD%ALP2*BMD%DT*BMD%G*BMD%HGU(I,J)*BMD%ETX(I,J)/BMD%DXU(I,J) * BMD%MSU(I,J)
          BMD%ETB(I-1,J) = BMD%ETB(I-1,J) - BMD%ALP2*BMD%DT*BMD%G*BMD%HGU(I,J)*BMD%ETX(I,J)/BMD%DXU(I,J) * BMD%MSU(I,J)
          BMD%ETX(I,J)   = 0.0
        ENDDO
       ENDDO

       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          BMD%UB(I+1,J) = BMD%UB(I+1,J) + BMD%DIV(I,J)/GRD%DX(I,J) * BMD%MST(I,J)
          BMD%UB(I  ,J) = BMD%UB(I  ,J) - BMD%DIV(I,J)/GRD%DX(I,J) * BMD%MST(I,J)
          BMD%VB(I,J+1) = BMD%VB(I,J+1) + BMD%DIV(I,J)/GRD%DY(I,J) * BMD%MST(I,J)
          BMD%VB(I,J  ) = BMD%VB(I,J  ) - BMD%DIV(I,J)/GRD%DY(I,J) * BMD%MST(I,J)
          BMD%DIV(I,J)  = 0.0
        ENDDO
       ENDDO

     ENDDO   ! KSTP


      DO J= 2,GRD%JM-1
       DO I= 2,GRD%IM-1
        BMD%BXBY(I,J) = BMD%BXBY(I,J) + BMD%RGH(I,J) * BMD%ALP1**2
        BMD%RGH(I,J)  = 0.0
        GRD%BX(I+1,J) = GRD%BX(I+1,J) - (BMD%DT**2)*BMD%BXBY(I,J)/GRD%DX(I,J) * BMD%MST(I,J)
        GRD%BX(I  ,J) = GRD%BX(I  ,J) + (BMD%DT**2)*BMD%BXBY(I,J)/GRD%DX(I,J) * BMD%MST(I,J)
        GRD%BY(I,J+1) = GRD%BY(I,J+1) - (BMD%DT**2)*BMD%BXBY(I,J)/GRD%DY(I,J) * BMD%MST(I,J)
        GRD%BY(I,J  ) = GRD%BY(I,J  ) + (BMD%DT**2)*BMD%BXBY(I,J)/GRD%DY(I,J) * BMD%MST(I,J)
        BMD%BXBY(I,J) = 0.0
       ENDDO
      ENDDO


END SUBROUTINE BAR_MOD_AD
