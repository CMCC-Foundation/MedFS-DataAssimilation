SUBROUTINE BAR_MOD

!-----------------------------------------------------------------------
!                                                                      !
! BAROTROPIC MODEL                                                     !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR

  IMPLICIT NONE

  INTEGER(I4)    :: I, J, KSTP
#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(J,I)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
      DO J= 2,GRD%JM-1
       DO I= 2,GRD%IM-1
        BMD%BXBY(I,J) = - (BMD%DT**2)*((GRD%BX(I+1,J)-GRD%BX(I,J))/GRD%DX(I,J) +               &
                                       (GRD%BY(I,J+1)-GRD%BY(I,J))/GRD%DY(I,J)) * BMD%MST(I,J)
         BMD%RGH(I,J) = BMD%BXBY(I,J) * BMD%ALP1**2
       ENDDO
      ENDDO
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif
       BMD%ETB(:,:) = 0.0  !BMD%ETM(:,:)
        BMD%UB(:,:) = 0.0  !BMD%UM(:,:)
        BMD%VB(:,:) = 0.0  !BMD%VM(:,:)
       BMD%ETN(:,:) = 0.0  !BMD%ETM(:,:)
        BMD%UN(:,:) = 0.0  !BMD%UM(:,:)
        BMD%VN(:,:) = 0.0  !BMD%VM(:,:)
       BMD%ETA(:,:) = 0.0  !BMD%ETM(:,:)
        BMD%UA(:,:) = 0.0  !BMD%UM(:,:)
        BMD%VA(:,:) = 0.0  !BMD%VM(:,:)

          BMD%ETX(:,:) = 0.0
          BMD%ETY(:,:) = 0.0
          BMD%DIV(:,:) = 0.0
           BMD%CU(:,:) = 0.0
           BMD%CV(:,:) = 0.0
          BMD%DUX(:,:) = 0.0
          BMD%DVX(:,:) = 0.0
          BMD%DUY(:,:) = 0.0
          BMD%DVY(:,:) = 0.0

         BMD%ETM(:,:) = 0.0
          BMD%UM(:,:) = 0.0
          BMD%VM(:,:) = 0.0


!---
! TIME LOOP
     DO KSTP = 1, BMD%NSTPS
#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(J,I)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
   DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          BMD%DIV(I,J) = ((BMD%UB(I+1,J)-BMD%UB(I,J))/GRD%DX(I,J) + (BMD%VB(I,J+1)-BMD%VB(I,J))/GRD%DY(I,J)) * BMD%MST(I,J)
        ENDDO
      ! ENDDO
      ! DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM
          BMD%ETX(I,J) = BMD%ALP2*BMD%DT*BMD%G*BMD%HGU(I,J)*(BMD%ETB(I,J)-BMD%ETB(I-1,J  ))/BMD%DXU(I,J) * BMD%MSU(I,J)
        ENDDO
       !ENDDO
       !DO J= 2,GRD%JM
        DO I= 2,GRD%IM-1
          BMD%ETY(I,J) = BMD%ALP2*BMD%DT*BMD%G*BMD%HGV(I,J)*(BMD%ETB(I,J)-BMD%ETB(I  ,J-1))/BMD%DYV(I,J) * BMD%MSV(I,J)
        ENDDO
       ENDDO
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif
       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
         BMD%RGH(I,J) = BMD%BXBY(I,J) - BMD%ETB(I,J) + BMD%DT*BMD%DIV(I,J)        &
                   - BMD%ALP1*BMD%DT*(BMD%ETX(I+1,J)-BMD%ETX(I,J))/GRD%DX(I,J)    &
                   - BMD%ALP1*BMD%DT*(BMD%ETY(I,J+1)-BMD%ETY(I,J))/GRD%DY(I,J)
        ENDDO
       ENDDO
       CALL INVRT( GRD%IM, GRD%JM, BMD%ETA, BMD%MST, BMD%RGH, BMD%A1, BMD%A2, BMD%A3, BMD%A4, BMD%A0,           &
                   BMD%BNM, BMD%OVR, BMD%RESEM, BMD%NCNT, BMD%ITR(KSTP) )

       DO J=1,GRD%JM-1
        DO I=2,GRD%IM-1
         BMD%CU(I,J) = -((BMD%VN(I,J  )*BMD%DXV(I,J  ) + BMD%VN(I-1,J  )*BMD%DXV(I-1,J  ))*GRD%F(I,J  ) +     &
                         (BMD%VN(I,J+1)*BMD%DXV(I,J+1) + BMD%VN(I-1,J+1)*BMD%DXV(I-1,J+1))*GRD%F(I,J+1))      &
                        * 0.25 / BMD%DXU(I,J)
        ENDDO
       ENDDO
       DO J=2,GRD%JM-1
        DO I=1,GRD%IM-1
         BMD%CV(I,J) =  ((BMD%UN(I  ,J)*BMD%DYU(I  ,J) + BMD%UN(I  ,J-1)*BMD%DYU(I  ,J-1))*GRD%F(I  ,J) +      &
                         (BMD%UN(I+1,J)*BMD%DYU(I+1,J) + BMD%UN(I+1,J-1)*BMD%DYU(I+1,J-1))*GRD%F(I+1,J))       &
                        * 0.25 / BMD%DYV(I,J)
        ENDDO
       ENDDO
          BMD%DUX(2:GRD%IM,:) = (BMD%UB(2:GRD%IM,:) - BMD%UB(1:GRD%IM-1,:))/BMD%DXU(2:GRD%IM,:)
          BMD%DVX(2:GRD%IM,:) = (BMD%VB(2:GRD%IM,:) - BMD%VB(1:GRD%IM-1,:))/BMD%DXU(2:GRD%IM,:)

          BMD%DUY(:,2:GRD%JM) = (BMD%UB(:,2:GRD%JM) - BMD%UB(:,1:GRD%JM-1))/BMD%DYV(:,2:GRD%JM)
          BMD%DVY(:,2:GRD%JM) = (BMD%VB(:,2:GRD%JM) - BMD%VB(:,1:GRD%JM-1))/BMD%DYV(:,2:GRD%JM)




       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          BMD%UA(I,J) = BMD%UB(I,J) - BMD%DT*(BMD%CU(I,J) +                                           &
                        BMD%ALP1*BMD%G*BMD%HGU(I,J)*(BMD%ETA(I,J)-BMD%ETA(I-1,J  ))/BMD%DXU(I,J) +    &
                        GRD%BX(I,J))*BMD%MSU(I,J)  -                                                  &
                        BMD%ETX(I,J) +                                                                &
                        BMD%DF1 * ((BMD%DUX(I+1,J)-BMD%DUX(I,J))/BMD%DXU(I,J) +                       &
                                   (BMD%DUY(I,J+1)-BMD%DUY(I,J))/BMD%DYU(I,J)) * BMD%MSU(I,J)
        ENDDO
       ENDDO
       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          BMD%VA(I,J) = BMD%VB(I,J) - BMD%DT*(BMD%CV(I,J) +                                           &
                        BMD%ALP1*BMD%G*BMD%HGV(I,J)*(BMD%ETA(I,J)-BMD%ETA(I  ,J-1))/BMD%DYV(I,J) +    &
                        GRD%BY(I,J))*BMD%MSV(I,J) -                                                   &
                        BMD%ETY(I,J) +                                                                &
                        BMD%DF1 * ((BMD%DVX(I+1,J)-BMD%DVX(I,J))/BMD%DXV(I,J) +                       &
                                   (BMD%DVY(I,J+1)-BMD%DVY(I,J))/BMD%DYV(I,J)) * BMD%MSV(I,J)
        ENDDO
       ENDDO

         DO J=1,GRD%JM
         DO I=2,GRD%IM
          BMD%DUX(I,J) = (BMD%UA(I,J) - BMD%UA(I-1,J))/BMD%DXU(I,J)
          BMD%DVX(I,J) = (BMD%VA(I,J) - BMD%VA(I-1,J))/BMD%DXV(I,J)
         ENDDO
         ENDDO
         DO J=2,GRD%JM
         DO I=1,GRD%IM
          BMD%DUY(I,J) = (BMD%UA(I,J) - BMD%UA(I,J-1))/BMD%DYU(I,J)
          BMD%DVY(I,J) = (BMD%VA(I,J) - BMD%VA(I,J-1))/BMD%DYV(I,J)
         ENDDO
         ENDDO

       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          BMD%UA(I,J) = BMD%UA(I,J) + BMD%DF2 * ((BMD%DUX(I+1,J)-BMD%DUX(I,J))/BMD%DXU(I,J) +      &
                                                 (BMD%DUY(I,J+1)-BMD%DUY(I,J))/BMD%DYU(I,J))*BMD%MSU(I,J)
        ENDDO
       ENDDO
       DO J= 2,GRD%JM-1
        DO I= 2,GRD%IM-1
          BMD%VA(I,J) = BMD%VA(I,J) + BMD%DF2 * ((BMD%DVX(I+1,J)-BMD%DVX(I,J))/BMD%DXV(I,J) +      &
                                                 (BMD%DVY(I,J+1)-BMD%DVY(I,J))/BMD%DYV(I,J))*BMD%MSV(I,J)
        ENDDO
       ENDDO

           BMD%UN(:,:) =  BMD%UN(:,:) + ( BMD%UB(:,:) + BMD%UA(:,:) - 2.0*BMD%UN(:,:) ) *0.05
           BMD%VN(:,:) =  BMD%VN(:,:) + ( BMD%VB(:,:) + BMD%VA(:,:) - 2.0*BMD%VN(:,:) ) *0.05

          BMD%ETB(:,:) = BMD%ETA(:,:)
           BMD%UB(:,:) =  BMD%UN(:,:)
           BMD%VB(:,:) =  BMD%VN(:,:)

           BMD%UN(:,:) =  BMD%UA(:,:)
           BMD%VN(:,:) =  BMD%VA(:,:)

          BMD%ETM(:,:) = BMD%ETM(:,:) + BMD%ETB(:,:)
           BMD%UM(:,:) =  BMD%UM(:,:) +  BMD%UB(:,:)
           BMD%VM(:,:) =  BMD%VM(:,:) +  BMD%VB(:,:)


       IF( MOD(KSTP,BMD%NSTPA).EQ.0 ) THEN

         BMD%ETM(:,:) = BMD%ETM(:,:)/BMD%NSTPA
          BMD%UM(:,:) =  BMD%UM(:,:)/BMD%NSTPA
          BMD%VM(:,:) =  BMD%VM(:,:)/BMD%NSTPA

         GRD%ETA(:,:) = BMD%ETM(:,:)

         BMD%ETM(:,:) = 0.0
          BMD%UM(:,:) = 0.0
          BMD%VM(:,:) = 0.0

       ENDIF


     ENDDO   ! KSTP
END SUBROUTINE BAR_MOD
