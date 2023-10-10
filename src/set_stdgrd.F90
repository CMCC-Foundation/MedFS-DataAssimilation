SUBROUTINE SET_STDGRD

!-----------------------------------------------------------------------
!                                                                      !
! DEFINE THE GRID                                                      !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR

  IMPLICIT NONE

  INTEGER(I4)    :: I, J, K
  REAL(R8)       :: PI, RADIUS

     GRD%IM = 800
     GRD%JM = 250
     GRD%KM = 48

     GRD%DLN = 56.0 / GRD%IM
     GRD%DLT = 17.0 / GRD%JM

     ALLOCATE ( GRD%REG(GRD%IM,GRD%JM))
     ALLOCATE ( GRD%MSK(GRD%IM,GRD%JM,GRD%KM))
     ALLOCATE ( GRD%TEM(GRD%IM,GRD%JM,GRD%KM), GRD%SAL(GRD%IM,GRD%JM,GRD%KM), GRD%ETA(GRD%IM,GRD%JM))
     ALLOCATE ( GRD%TEM_AD(GRD%IM,GRD%JM,GRD%KM), GRD%SAL_AD(GRD%IM,GRD%JM,GRD%KM), GRD%ETA_AD(GRD%IM,GRD%JM))
     ALLOCATE ( GRD%LON(GRD%IM,GRD%JM), GRD%LAT(GRD%IM,GRD%JM), GRD%DEP(GRD%KM))
     ALLOCATE ( GRD%DX(GRD%IM,GRD%JM), GRD%DY(GRD%IM,GRD%JM), GRD%DXDY(GRD%IM,GRD%JM))


     DO J=1,GRD%JM
     DO I=1,GRD%IM
      GRD%LON(I,J) =  -10. + (I-1)*GRD%DLN
      GRD%LAT(I,J) = 28. + (J-1)*GRD%DLT
     ENDDO
     ENDDO

     GRD%MSK(:,:,:) = 1.0

     DO J=1,GRD%JM
     DO I=1,GRD%IM
      GRD%DY(I,J) = 5.E3
      GRD%DX(I,J) = 5.E3
      GRD%DXDY(I,J) = GRD%DY(I,J) * GRD%DX(I,J)
     ENDDO
     ENDDO
      GRD%ADXDY = SUM(GRD%DXDY) / (GRD%IM*GRD%JM)

      GRD%REG(:,:) = 1

     DO K=1,GRD%KM
      GRD%DEP(K) =  0. + (K-1)*10.
     ENDDO
      GRD%NPS = GRD%IM*GRD%JM*GRD%KM

RETURN
END SUBROUTINE SET_STDGRD
