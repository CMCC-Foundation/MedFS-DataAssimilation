SUBROUTINE GET_VEL_AD

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE VERTICAL INTEGRAL OF BOUYANCY GRADIENT                     !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR

  IMPLICIT NONE

  REAL(R8), DIMENSION (GRD%IM,GRD%JM)  :: UD, VD

  INTEGER(I4)    :: K

      UD(:,:) = 0.0
      VD(:,:) = 0.0

     DO K=GRD%KM,1,-1


      VD(1:GRD%IM-1,2:GRD%JM)   = VD(1:GRD%IM-1,2:GRD%JM)   + GRD%VVL_AD(1:GRD%IM-1,2:GRD%JM,K)*0.25    &
                                * GRD%MSK(1:GRD%IM-1,2:GRD%JM,K)*GRD%MSK(1:GRD%IM-1,1:GRD%JM-1,K)
      VD(1:GRD%IM-1,1:GRD%JM-1) = VD(1:GRD%IM-1,1:GRD%JM-1) + GRD%VVL_AD(1:GRD%IM-1,2:GRD%JM,K)*0.25    &
                                * GRD%MSK(1:GRD%IM-1,2:GRD%JM,K)*GRD%MSK(1:GRD%IM-1,1:GRD%JM-1,K)
      VD(2:GRD%IM,2:GRD%JM)     = VD(2:GRD%IM,2:GRD%JM)     + GRD%VVL_AD(1:GRD%IM-1,2:GRD%JM,K)*0.25    &
                                * GRD%MSK(2:GRD%IM,2:GRD%JM,K)*GRD%MSK(2:GRD%IM,1:GRD%JM-1,K)
      VD(2:GRD%IM,1:GRD%JM-1)   = VD(2:GRD%IM,1:GRD%JM-1)   + GRD%VVL_AD(1:GRD%IM-1,2:GRD%JM,K)*0.25    &
                                * GRD%MSK(2:GRD%IM,2:GRD%JM,K)*GRD%MSK(2:GRD%IM,1:GRD%JM-1,K)

      GRD%VVL_AD(1:GRD%IM-1,2:GRD%JM,K) = 0.0

      UD(2:GRD%IM,1:GRD%JM-1)   = UD(2:GRD%IM,1:GRD%JM-1)   + GRD%UVL_AD(2:GRD%IM,1:GRD%JM-1,K)*0.25    &
                                * GRD%MSK(2:GRD%IM,1:GRD%JM-1,K)*GRD%MSK(1:GRD%IM-1,1:GRD%JM-1,K)
      UD(1:GRD%IM-1,1:GRD%JM-1) = UD(1:GRD%IM-1,1:GRD%JM-1) + GRD%UVL_AD(2:GRD%IM,1:GRD%JM-1,K)*0.25    &
                                * GRD%MSK(2:GRD%IM,1:GRD%JM-1,K)*GRD%MSK(1:GRD%IM-1,1:GRD%JM-1,K)
      UD(2:GRD%IM,2:GRD%JM)     = UD(2:GRD%IM,2:GRD%JM)     + GRD%UVL_AD(2:GRD%IM,1:GRD%JM-1,K)*0.25    &
                                * GRD%MSK(2:GRD%IM,2:GRD%JM,K)*GRD%MSK(1:GRD%IM-1,2:GRD%JM,K)
      UD(1:GRD%IM-1,2:GRD%JM)   = UD(1:GRD%IM-1,2:GRD%JM)   + GRD%UVL_AD(2:GRD%IM,1:GRD%JM-1,K)*0.25    &
                                * GRD%MSK(2:GRD%IM,2:GRD%JM,K)*GRD%MSK(1:GRD%IM-1,2:GRD%JM,K)

      GRD%UVL_AD(2:GRD%IM,1:GRD%JM-1,K) = 0.0

      GRD%B_Y(:,2:GRD%JM,K)    = GRD%B_Y(:,2:GRD%JM,K)    - UD(:,2:GRD%JM)       &
                              / GRD%DY(:,2:GRD%JM) * GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K) / GRD%F(:,2:GRD%JM)
      GRD%ETA_AD(:,1:GRD%JM-1) = GRD%ETA_AD(:,1:GRD%JM-1) + UD(:,2:GRD%JM)*9.81       &
                              / GRD%DY(:,2:GRD%JM) * GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K) / GRD%F(:,2:GRD%JM)
      GRD%ETA_AD(:,2:GRD%JM)   = GRD%ETA_AD(:,2:GRD%JM)   - UD(:,2:GRD%JM)*9.81       &
                              / GRD%DY(:,2:GRD%JM) * GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K) / GRD%F(:,2:GRD%JM)

      GRD%B_X(2:GRD%IM,:,K)    = GRD%B_X(2:GRD%IM,:,K)    + VD(2:GRD%IM,:)       &
                              / GRD%DX(2:GRD%IM,:) * GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K) / GRD%F(2:GRD%IM,:)
      GRD%ETA_AD(1:GRD%IM-1,:) = GRD%ETA_AD(1:GRD%IM-1,:) - VD(2:GRD%IM,:)*9.81       &
                              / GRD%DX(2:GRD%IM,:) * GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K) / GRD%F(2:GRD%IM,:)
      GRD%ETA_AD(2:GRD%IM,:)   = GRD%ETA_AD(2:GRD%IM,:)   + VD(2:GRD%IM,:)*9.81       &
                              / GRD%DX(2:GRD%IM,:) * GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K) / GRD%F(2:GRD%IM,:)

      UD(:,:) = 0.0
      VD(:,:) = 0.0

     ENDDO

     GRD%UVL_AD(:,:,:) = 0.0
     GRD%VVL_AD(:,:,:) = 0.0

      GRD%B_Y(:,2:GRD%JM,1) = GRD%B_Y(:,2:GRD%JM,1) * 0.5
      GRD%B_X(2:GRD%IM,:,1) = GRD%B_X(2:GRD%IM,:,1) * 0.5

     DO K=2,GRD%KM
      GRD%B_Y(:,2:GRD%JM,K-1) = GRD%B_Y(:,2:GRD%JM,K-1) + GRD%B_Y(:,2:GRD%JM,K) * 0.5
      GRD%B_Y(:,2:GRD%JM,K)   = GRD%B_Y(:,2:GRD%JM,K) * 0.5
      GRD%B_X(2:GRD%IM,:,K-1) = GRD%B_X(2:GRD%IM,:,K-1) + GRD%B_X(2:GRD%IM,:,K) * 0.5
      GRD%B_X(2:GRD%IM,:,K)   = GRD%B_X(2:GRD%IM,:,K) * 0.5
     ENDDO

END SUBROUTINE GET_VEL_AD
