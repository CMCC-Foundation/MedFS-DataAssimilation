SUBROUTINE DIV_DMP_AD

!-----------------------------------------------------------------------
!                                                                      !
! DIVERGENCE DAMPING (ADJOINT)
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR
  USE MYFRTPROF, ONLY : MYFRTPROF_WALL

  IMPLICIT NONE

  REAL(R8), DIMENSION (GRD%IM,GRD%JM)  :: DIV

  INTEGER(I4)    :: K, KDIV

  CALL MYFRTPROF_WALL('DIV_DMP_AD: ADJOINT OF DIVERGENCE DAMPING',0)

      DIV(:,:) = 0.0

     DO K=GRD%KM,1,-1

      DO KDIV = 1,100

       DIV(1:GRD%IM-1,2:GRD%JM-1) = DIV(1:GRD%IM-1,2:GRD%JM-1) + GRD%VVL(1:GRD%IM-1,2:GRD%JM-1,K)*0.2 * GRD%ADXDY**2    &
                                  / GRD%DY(1:GRD%IM-1,2:GRD%JM-1)                                                       &
                                  * GRD%MSK(1:GRD%IM-1,2:GRD%JM-1,K)*GRD%MSK(1:GRD%IM-1,1:GRD%JM-2,K)
       DIV(1:GRD%IM-1,1:GRD%JM-2) = DIV(1:GRD%IM-1,1:GRD%JM-2) - GRD%VVL(1:GRD%IM-1,2:GRD%JM-1,K)*0.2 * GRD%ADXDY**2    &
                                  / GRD%DY(1:GRD%IM-1,2:GRD%JM-1)                                                       &
                                  * GRD%MSK(1:GRD%IM-1,2:GRD%JM-1,K)*GRD%MSK(1:GRD%IM-1,1:GRD%JM-2,K)
       DIV(2:GRD%IM-1,1:GRD%JM-1) = DIV(2:GRD%IM-1,1:GRD%JM-1) + GRD%UVL(2:GRD%IM-1,1:GRD%JM-1,K)*0.2 * GRD%ADXDY**2    &
                                  / GRD%DX(2:GRD%IM-1,1:GRD%JM-1)                                                       &
                                  * GRD%MSK(2:GRD%IM-1,1:GRD%JM-1,K)*GRD%MSK(1:GRD%IM-2,1:GRD%JM-1,K)
       DIV(1:GRD%IM-2,1:GRD%JM-1) = DIV(1:GRD%IM-2,1:GRD%JM-1) - GRD%UVL(2:GRD%IM-1,1:GRD%JM-1,K)*0.2 * GRD%ADXDY**2    &
                                  / GRD%DX(2:GRD%IM-1,1:GRD%JM-1)                                                       &
                                  * GRD%MSK(2:GRD%IM-1,1:GRD%JM-1,K)*GRD%MSK(1:GRD%IM-2,1:GRD%JM-1,K)

       GRD%UVL(2:GRD%IM  ,1:GRD%JM-1,K) = GRD%UVL(2:GRD%IM  ,1:GRD%JM-1,K) + DIV(1:GRD%IM-1,1:GRD%JM-1)    &
                                          * GRD%DY(2:GRD%IM  ,1:GRD%JM-1)                                  &
                                          / GRD%DX(1:GRD%IM-1,1:GRD%JM-1) / GRD%DY(1:GRD%IM-1,1:GRD%JM-1)
       GRD%UVL(1:GRD%IM-1,1:GRD%JM-1,K) = GRD%UVL(1:GRD%IM-1,1:GRD%JM-1,K) - DIV(1:GRD%IM-1,1:GRD%JM-1)    &
                                          * GRD%DY(1:GRD%IM-1,1:GRD%JM-1)                                  &
                                          / GRD%DX(1:GRD%IM-1,1:GRD%JM-1) / GRD%DY(1:GRD%IM-1,1:GRD%JM-1)
       GRD%VVL(1:GRD%IM-1,2:GRD%JM  ,K) = GRD%VVL(1:GRD%IM-1,2:GRD%JM  ,K) + DIV(1:GRD%IM-1,1:GRD%JM-1)    &
                                          * GRD%DX(1:GRD%IM-1,2:GRD%JM  )                                  &
                                          / GRD%DX(1:GRD%IM-1,1:GRD%JM-1) / GRD%DY(1:GRD%IM-1,1:GRD%JM-1)
       GRD%VVL(1:GRD%IM-1,1:GRD%JM-1,K) = GRD%VVL(1:GRD%IM-1,1:GRD%JM-1,K) - DIV(1:GRD%IM-1,1:GRD%JM-1)    &
                                          * GRD%DX(1:GRD%IM-1,1:GRD%JM-1)                                  &
                                          / GRD%DX(1:GRD%IM-1,1:GRD%JM-1) / GRD%DY(1:GRD%IM-1,1:GRD%JM-1)
       DIV(:,:) = 0.0

      ENDDO

     ENDDO

  CALL MYFRTPROF_WALL('DIV_DMP_AD: ADJOINT OF DIVERGENCE DAMPING',1)
END SUBROUTINE DIV_DMP_AD
