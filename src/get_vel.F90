SUBROUTINE GET_VEL

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
  REAL(R8), DIMENSION (GRD%IM,GRD%JM)  :: DIV

  INTEGER(I4)    :: K, KDIV

     DO K=GRD%KM,2,-1
      GRD%B_X(2:GRD%IM,:,K) = ( GRD%B_X(2:GRD%IM,:,K) + GRD%B_X(2:GRD%IM,:,K-1) ) * 0.5
      GRD%B_Y(:,2:GRD%JM,K) = ( GRD%B_Y(:,2:GRD%JM,K) + GRD%B_Y(:,2:GRD%JM,K-1) ) * 0.5
     ENDDO
      GRD%B_X(2:GRD%IM,:,1) = GRD%B_X(2:GRD%IM,:,1) * 0.5
      GRD%B_Y(:,2:GRD%JM,1) = GRD%B_Y(:,2:GRD%JM,1) * 0.5

     GRD%UVL(:,:,:) = 0.0
     GRD%VVL(:,:,:) = 0.0

     DO K=1,GRD%KM

      UD(:,:) = 0.0
      VD(:,:) = 0.0

      VD(2:GRD%IM,:) =   ( (GRD%ETA(2:GRD%IM,:)-GRD%ETA(1:GRD%IM-1,:))*9.81 + GRD%B_X(2:GRD%IM,:,K) )               &
                           / GRD%DX(2:GRD%IM,:) * GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K) / GRD%F(2:GRD%IM,:)
      UD(:,2:GRD%JM) = - ( (GRD%ETA(:,2:GRD%JM)-GRD%ETA(:,1:GRD%JM-1))*9.81 + GRD%B_Y(:,2:GRD%JM,K) )               &
                           / GRD%DY(:,2:GRD%JM) * GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K) / GRD%F(:,2:GRD%JM)
      GRD%UVL(2:GRD%IM,1:GRD%JM-1,K) =                     &
                         ( UD(2:GRD%IM,1:GRD%JM-1)+UD(1:GRD%IM-1,1:GRD%JM-1)+UD(2:GRD%IM,2:GRD%JM)+UD(1:GRD%IM-1,2:GRD%JM) )*0.25  &
                         * GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K)
      GRD%VVL(1:GRD%IM-1,2:GRD%JM,K) =                     &
                         ( VD(1:GRD%IM-1,2:GRD%JM)+VD(1:GRD%IM-1,1:GRD%JM-1)+VD(2:GRD%IM,2:GRD%JM)+VD(2:GRD%IM,1:GRD%JM-1) )*0.25  &
                         * GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K)

     ENDDO


END SUBROUTINE GET_VEL
