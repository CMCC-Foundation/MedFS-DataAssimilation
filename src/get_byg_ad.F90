SUBROUTINE GET_BYG_AD

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE VERTICAL INTEGRAL OF BOUYANCY GRADIENT (ADJOINT)           !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR
  USE MOD_BMD
  USE OCEANTOOLS
  USE MYFRTPROF, ONLY : MYFRTPROF_WALL

  IMPLICIT NONE

  INTEGER(I4)    :: K,J,I
  REAL(R8)       :: SAD, TAD
CALL MYFRTPROF_WALL('GET_BYG_AD: ADJOINT OF BUOYANCY GRADIENT',0)

      GRD%DNS(:,:,:) = 0.0

      GRD%BX(2:GRD%IM,:) = GRD%BX(2:GRD%IM,:)/GRD%DX(2:GRD%IM,:)
      GRD%BY(:,2:GRD%JM) = GRD%BY(:,2:GRD%JM)/GRD%DY(:,2:GRD%JM)

     DO K=1,GRD%KM
      GRD%B_X(2:GRD%IM,:,K) = GRD%B_X(2:GRD%IM,:,K) + GRD%BX(2:GRD%IM,:)*GRD%DZ(K) * GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K)
      GRD%B_Y(:,2:GRD%JM,K) = GRD%B_Y(:,2:GRD%JM,K) + GRD%BY(:,2:GRD%JM)*GRD%DZ(K) * GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K)
     ENDDO

     DO K=GRD%KM,2,-1
      GRD%DNS(:,2:GRD%JM  ,K) = GRD%DNS(:,2:GRD%JM  ,K) +               &
                                GRD%B_Y(:,2:GRD%JM,K)*GRD%DZ(K)*GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K)
      GRD%DNS(:,1:GRD%JM-1,K) = GRD%DNS(:,1:GRD%JM-1,K) -               &
                                GRD%B_Y(:,2:GRD%JM,K)*GRD%DZ(K)*GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K)
      GRD%B_Y(:,2:GRD%JM,K-1) = GRD%B_Y(:,2:GRD%JM,K-1) + GRD%B_Y(:,2:GRD%JM,K)
      GRD%B_Y(:,2:GRD%JM,K)   = 0.0
     ENDDO
      GRD%DNS(:,2:GRD%JM  ,1) = GRD%DNS(:,2:GRD%JM  ,1) +               &
                                GRD%B_Y(:,2:GRD%JM,1)*GRD%DZ(1)*GRD%MSK(:,2:GRD%JM,1)*GRD%MSK(:,1:GRD%JM-1,1)
      GRD%DNS(:,1:GRD%JM-1,1) = GRD%DNS(:,1:GRD%JM-1,1) -               &
                                GRD%B_Y(:,2:GRD%JM,1)*GRD%DZ(1)*GRD%MSK(:,2:GRD%JM,1)*GRD%MSK(:,1:GRD%JM-1,1)

     DO K=GRD%KM,2,-1
      GRD%DNS(2:GRD%IM  ,:,K) = GRD%DNS(2:GRD%IM  ,:,K) +               &
                                GRD%B_X(2:GRD%IM,:,K)*GRD%DZ(K)*GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K)
      GRD%DNS(1:GRD%IM-1,:,K) = GRD%DNS(1:GRD%IM-1,:,K) -               &
                                GRD%B_X(2:GRD%IM,:,K)*GRD%DZ(K)*GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K)
      GRD%B_X(2:GRD%IM,:,K-1) = GRD%B_X(2:GRD%IM,:,K-1) + GRD%B_X(2:GRD%IM,:,K)
      GRD%B_X(2:GRD%IM,:,K)   = 0.0
     ENDDO
      GRD%DNS(2:GRD%IM  ,:,1) = GRD%DNS(2:GRD%IM  ,:,1) +               &
                                GRD%B_X(2:GRD%IM,:,1)*GRD%DZ(K)*GRD%MSK(2:GRD%IM,:,1)*GRD%MSK(1:GRD%IM-1,:,1)
      GRD%DNS(1:GRD%IM-1,:,1) = GRD%DNS(1:GRD%IM-1,:,1) -               &
                                GRD%B_X(2:GRD%IM,:,1)*GRD%DZ(K)*GRD%MSK(2:GRD%IM,:,1)*GRD%MSK(1:GRD%IM-1,:,1)

       GRD%DNS(:,:,:) = GRD%DNS(:,:,:) *9.81/1020._R8

!      GRD%TEM_AD(:,:,:) = GRD%TEM_AD(:,:,:) - 0.24*GRD%DNS(:,:,:) * 9.81/1025. * GRD%MSR(:,:,:)
!      GRD%SAL_AD(:,:,:) = GRD%SAL_AD(:,:,:) + 0.74*GRD%DNS(:,:,:) * 9.81/1025. * GRD%MSR(:,:,:)

 IF( BMD_NNEOS .EQ. 1 ) THEN
     DO K=1,GRD%KM
      GRD%TEM_AD(:,:,K) = GRD%TEM_AD(:,:,K) - 0.24*GRD%DNS(:,:,K) * GRD%MSK(:,:,K)
      GRD%SAL_AD(:,:,K) = GRD%SAL_AD(:,:,K) + 0.74*GRD%DNS(:,:,K) * GRD%MSK(:,:,K)
     ENDDO
   ELSEIF( BMD_NNEOS .EQ. 2 ) THEN
     DO K=1,GRD%KM
        DO J=1,GRD%JM
           DO I=1,GRD%IM
              CALL RHO_UNESCOAD(GRD%DNS(I,J,K),BMD_SALB(I,J,K),BMD_TEMB(I,J,K),SAD, TAD)
              GRD%TEM_AD(I,J,K) = GRD%TEM_AD(I,J,K) + TAD*GRD%MSK(I,J,K)
              GRD%SAL_AD(I,J,K) = GRD%SAL_AD(I,J,K) + SAD*GRD%MSK(I,J,K)
           ENDDO
        ENDDO
     ENDDO
ELSE
       CALL ABOR1('GET_BYD_AD: UNSUPPORTED E.O.S. OPTION')
   ENDIF




CALL MYFRTPROF_WALL('GET_BYG_AD: ADJOINT OF BUOYANCY GRADIENT',1)
END SUBROUTINE GET_BYG_AD
