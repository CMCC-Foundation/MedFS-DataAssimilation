SUBROUTINE INI_BMD

!-----------------------------------------------------------------------
!                                                                      !
! INITIALISE THE BAROTROPIC MODEL                                      !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR
  USE MOD_BMD
  USE MYNETCDF 
  USE IOUNITS 
  IMPLICIT NONE

  INTEGER(I4)    :: I, J

     BMD%G     = 9.80665_R8
     BMD%NSTP  = INT( 24._R8 * 3600._R8 / BMD%DT )
     BMD%NSTPS = INT( BMD%NDY * BMD%NSTP )
     BMD%NSTPA = INT( BMD%ADY * BMD%NSTP )
     BMD%ALP2  = 1.0_R8 - BMD%ALP1

     BMD%DF1 = BMD%FC1 * GRD%ADXDY**2
     BMD%DF2 = BMD%FC2 * GRD%ADXDY**2

     ALLOCATE ( BMD%ITR(BMD%NSTPS) )
     ALLOCATE ( BMD%MST(GRD%IM,GRD%JM), BMD%MSU(GRD%IM,GRD%JM), BMD%MSV(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%HGT(GRD%IM,GRD%JM), BMD%HGU(GRD%IM,GRD%JM), BMD%HGV(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%DXU(GRD%IM,GRD%JM), BMD%DXV(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%DYU(GRD%IM,GRD%JM), BMD%DYV(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%A1(GRD%IM,GRD%JM), BMD%A2(GRD%IM,GRD%JM), BMD%A3(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%A4(GRD%IM,GRD%JM), BMD%A0(GRD%IM,GRD%JM), BMD%A00(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%BX(GRD%IM,GRD%JM), BMD%BY(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%B_X(GRD%IM,GRD%JM,GRD%KM), BMD%B_Y(GRD%IM,GRD%JM,GRD%KM))
     ALLOCATE ( BMD%DNS(GRD%IM,GRD%JM,GRD%KM) )
     ALLOCATE ( BMD%BXBY(GRD%IM,GRD%JM), BMD%RGH(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%ETB(GRD%IM,GRD%JM), BMD%UB(GRD%IM,GRD%JM), BMD%VB(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%ETN(GRD%IM,GRD%JM), BMD%UN(GRD%IM,GRD%JM), BMD%VN(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%ETA(GRD%IM,GRD%JM), BMD%UA(GRD%IM,GRD%JM), BMD%VA(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%ETM(GRD%IM,GRD%JM), BMD%UM(GRD%IM,GRD%JM), BMD%VM(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%DIV(GRD%IM,GRD%JM), BMD%CU(GRD%IM,GRD%JM), BMD%CV(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%DUX(GRD%IM,GRD%JM), BMD%DUY(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%DVX(GRD%IM,GRD%JM), BMD%DVY(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%ETX(GRD%IM,GRD%JM), BMD%ETY(GRD%IM,GRD%JM))
     ALLOCATE ( BMD%BFU(GRD%IM,GRD%JM), BMD%BFV(GRD%IM,GRD%JM))

     bmd%mst(:,:) = 0._R8
     bmd%msu(:,:) = 0._R8
     bmd%msv(:,:) = 0._R8
     bmd%hgt(:,:) = 0._R8
     bmd%hgu(:,:) = 0._R8
     bmd%hgv(:,:) = 0._R8
     bmd%dxu(:,:) = 0._R8
     bmd%dyv(:,:) = 0._R8
     bmd%dyu(:,:) = 0._R8
     bmd%dxv(:,:) = 0._R8
     bmd%a1 (:,:) = 0._R8
     bmd%a2 (:,:) = 0._R8
     bmd%a3 (:,:) = 0._R8
     bmd%a4 (:,:) = 0._R8
     bmd%a0 (:,:) = 0._R8
     bmd%a00(:,:) = 0._R8

     BMD%ITR(:) = 0

     DO J=1,GRD%JM
       DO I=1,GRD%IM
         IF(GRD%HGT(I,J).GT.0._R8) THEN
           BMD%HGT(I,J) = GRD%HGT(I,J)
         ELSE
           BMD%HGT(I,J) = 0._R8
         ENDIF
       ENDDO
     ENDDO

     BMD%HGT(1,:) = 0._R8
     BMD%HGT(:,1) = 0._R8
     BMD%HGT(GRD%IM,:) = 0._R8
     BMD%HGT(:,GRD%JM) = 0._R8

!-------------------------------------------------

     DO J=1,GRD%JM
       DO I=1,GRD%IM
         IF(BMD%HGT(I,J).GT.0._R8) THEN
           BMD%MST(I,J) = 1._R8
         ELSE
           BMD%MST(I,J) = 0._R8
         ENDIF
       ENDDO
     ENDDO

     BMD%BNM = 0._R8
     DO J=2,GRD%JM-1
       DO I=2,GRD%IM-1
         IF(BMD%MST(I,J).EQ.1.) BMD%BNM = BMD%BNM + 1.
       ENDDO
     ENDDO

     DO J=1,GRD%JM
       DO I=2,GRD%IM
         BMD%HGU(I,J) = MIN(BMD%HGT(I,J),BMD%HGT(I-1,J))
         BMD%DXU(I,J) = (GRD%DX(I,J)+GRD%DX(I-1,J))*0.5_R8
         BMD%DYU(I,J) = (GRD%DY(I,J)+GRD%DY(I-1,J))*0.5_R8
         BMD%MSU(I,J) =  BMD%MST(I,J)*BMD%MST(I-1,J)
       ENDDO
     ENDDO
     BMD%DXU(1,:) = BMD%DXU(2,:)
     BMD%DYU(1,:) = BMD%DYU(2,:)
     BMD%MSU(1,:) = 0._R8

     DO J=2,GRD%JM
       DO I=1,GRD%IM
         BMD%HGV(I,J) = MIN(BMD%HGT(I,J),BMD%HGT(I,J-1))
         BMD%DXV(I,J) = (GRD%DX(I,J)+GRD%DX(I,J-1))*0.5_R8
         BMD%DYV(I,J) = (GRD%DY(I,J)+GRD%DY(I,J-1))*0.5_R8
         BMD%MSV(I,J) =  BMD%MST(I,J)*BMD%MST(I,J-1)
       ENDDO
     ENDDO
     BMD%DXV(:,1) = BMD%DXV(:,2)
     BMD%DYV(:,1) = BMD%DYV(:,2)
     BMD%MSV(:,1) = 0._R8

     DO J= 2,GRD%JM-1
       DO I= 2,GRD%IM-1
         BMD%A1(I,J) = BMD%ALP1**2 *(BMD%DT**2)*BMD%G*BMD%HGU(I+1,J)/GRD%DX(I,J)**2*BMD%MSU(I+1,J)
         BMD%A2(I,J) = BMD%ALP1**2 *(BMD%DT**2)*BMD%G*BMD%HGU(I  ,J)/GRD%DX(I,J)**2*BMD%MSU(I  ,J)
         BMD%A3(I,J) = BMD%ALP1**2 *(BMD%DT**2)*BMD%G*BMD%HGV(I,J+1)/GRD%DY(I,J)**2*BMD%MSV(I,J+1)
         BMD%A4(I,J) = BMD%ALP1**2 *(BMD%DT**2)*BMD%G*BMD%HGV(I,J  )/GRD%DY(I,J)**2*BMD%MSV(I,J  )
         BMD%A0(I,J) = (BMD%A1(I,J)+BMD%A2(I,J)+BMD%A3(I,J)+BMD%A4(I,J) + 1._R8)
        BMD%A00(I,J) = (BMD%A1(I,J)+BMD%A2(I,J)+BMD%A3(I,J)+BMD%A4(I,J)) *BMD%MST(I,J)
       ENDDO
     ENDDO

     ! BOTTOM FRICTION

     DO J= 1,GRD%JM
       DO I= 2,GRD%IM
         bmd%bfu(i,j) = 0.1_R8/max(0.1_r8,bmd%hgu(i,j)) * bmd%msu(i,j)
         bmd%bfu(i,j) = 0.01_R8* bmd%msu(i,j)
       ENDDO
     ENDDO

     DO J= 2,GRD%JM
       DO I= 1,GRD%IM
         bmd%bfv(i,j) = 0.1_R8/max(0.1_r8,bmd%hgv(i,j)) * bmd%msv(i,j)
         bmd%bfv(i,j) = 0.01_R8* bmd%msv(i,j)
       ENDDO
     ENDDO
   IF( BMD_NNEOS .EQ. 2 ) THEN
    ALLOCATE( BMD_TEMB(GRD%IM,GRD%JM,GRD%KM) )
    ALLOCATE( BMD_SALB(GRD%IM,GRD%JM,GRD%KM) )
    CALL GETNCVAR('BACKGROUND.nc','votemper',GRD%IM,GRD%JM,GRD%KM,BMD_TEMB)
    CALL GETNCVAR('BACKGROUND.nc','vosaline',GRD%IM,GRD%JM,GRD%KM,BMD_SALB)
  ENDIF

END SUBROUTINE INI_BMD
