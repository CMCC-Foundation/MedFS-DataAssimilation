SUBROUTINE DIV_DMP

!-----------------------------------------------------------------------
!                                                                      !
! DIVERGENCE DAMPING
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR
  USE MYFRTPROF, ONLY :MYFRTPROF_WALL

  IMPLICIT NONE

  REAL(R8), DIMENSION (GRD%IM,GRD%JM)  :: DIV

  INTEGER(I4)    :: K, KDIV

  CALL MYFRTPROF_WALL('DIV_DMP: DIVERGENCE DAMPING',0)


#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(K,KDIV,DIV)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
     DO K=1,GRD%KM

      DO KDIV = 1,100
       DIV(:,:) = 0.0
       DIV(1:GRD%IM-1,1:GRD%JM-1) = ( GRD%UVL(2:GRD%IM  ,1:GRD%JM-1,K)*&
                                  & GRD%DY(2:GRD%IM  ,1:GRD%JM-1) -       &
                                      GRD%UVL(1:GRD%IM-1,1:GRD%JM-1,K)*GRD%DY(1:GRD%IM-1,1:GRD%JM-1) +       &
                                      GRD%VVL(1:GRD%IM-1,2:GRD%JM  ,K)*GRD%DX(1:GRD%IM-1,2:GRD%JM  ) -       &
                                      GRD%VVL(1:GRD%IM-1,1:GRD%JM-1,K)*GRD%DX(1:GRD%IM-1,1:GRD%JM-1) )       &
                                    / GRD%DX(1:GRD%IM-1,1:GRD%JM-1) / GRD%DY(1:GRD%IM-1,1:GRD%JM-1)
       GRD%UVL(2:GRD%IM-1,1:GRD%JM-1,K) = GRD%UVL(2:GRD%IM-1,1:GRD%JM-1,K)  +                                &
                           0.2  * GRD%ADXDY**2 * (DIV(2:GRD%IM-1,1:GRD%JM-1) - DIV(1:GRD%IM-2,1:GRD%JM-1))    &
                           / GRD%DX(2:GRD%IM-1,1:GRD%JM-1) * GRD%MSK(2:GRD%IM-1,1:GRD%JM-1,K)*GRD%MSK(1:GRD%IM-2,1:GRD%JM-1,K)
     GRD%VVL(1:GRD%IM-1,2:GRD%JM-1,K) = GRD%VVL(1:GRD%IM-1,2:GRD%JM-1,K)  +                                &
                           0.2  * GRD%ADXDY**2 * (DIV(1:GRD%IM-1,2:GRD%JM-1) - DIV(1:GRD%IM-1,1:GRD%JM-2))    &
                           / GRD%DY(1:GRD%IM-1,2:GRD%JM-1) * GRD%MSK(1:GRD%IM-1,2:GRD%JM-1,K)*GRD%MSK(1:GRD%IM-1,1:GRD%JM-2,K)
      ENDDO

     ENDDO
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif
  CALL MYFRTPROF_WALL('DIV_DMP: DIVERGENCE DAMPING',1)
END SUBROUTINE DIV_DMP
