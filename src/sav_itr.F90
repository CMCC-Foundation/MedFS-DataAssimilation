SUBROUTINE SAV_ITR

!-----------------------------------------------------------------------
!                                                                      !
! SAVE THE RESULT ON THE COARSE GRID                                   !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE DRV_STR
 USE OBS_STR
 USE GRD_STR
 USE EOF_STR
 USE CTL_STR
 USE BMD_STR

 IMPLICIT NONE

 INTEGER(I4)    :: I,J,K, KK
 REAL(R8)    :: TMN,SMN, SPD,ONEN

! ---
! SAVE GRID DIMENSIONS

   DRV%IM = GRD%IM
   DRV%JM = GRD%JM

   ALLOCATE ( DRV%RO(DRV%IM,DRV%JM,ROS%NEOF), DRV%RO_AD(DRV%IM,DRV%JM,ROS%NEOF))
   ALLOCATE ( DRV%MSK(DRV%IM,DRV%JM))

! ---
! SAVE EIGENVALUES

   DRV%RO   (:,:,:) = GRD%RO   (:,:,:)
   DRV%RO_AD(:,:,:) = GRD%RO_AD(:,:,:)
   DRV%MSK  (:,:)   = GRD%MSR  (:,:,1)

! ---
! DEALLOCATE EVERITHING RELATED TO THE OLD GRID

! GRID STRUCTURE
     DEALLOCATE ( GRD%REG)
     DEALLOCATE ( GRD%MSK)
     DEALLOCATE ( GRD%HGT)
     DEALLOCATE ( GRD%F)
     DEALLOCATE ( GRD%TEM, GRD%SAL)
     DEALLOCATE ( GRD%UVL, GRD%VVL)
     DEALLOCATE ( GRD%UVL_AD, GRD%VVL_AD)
     DEALLOCATE ( GRD%B_X, GRD%B_Y)
     DEALLOCATE ( GRD%DNS)
     DEALLOCATE ( GRD%BX, GRD%BY)
     DEALLOCATE ( GRD%ETA)
     DEALLOCATE ( GRD%MDT, GRD%SLA)
     DEALLOCATE ( GRD%TEM_AD, GRD%SAL_AD)
     DEALLOCATE ( GRD%ETA_AD)
     DEALLOCATE ( GRD%LON, GRD%LAT, GRD%DEP)
     DEALLOCATE ( GRD%DX, GRD%DY, GRD%DZ)
     DEALLOCATE ( GRD%DXDY)
     DEALLOCATE ( GRD%MSR )
! OBSERVATIONAL VECTOR
     DEALLOCATE ( OBS%INC, OBS%AMO, OBS%RES)
     DEALLOCATE ( OBS%ERR, OBS%GRA)
! COVARIANCES STRUCTURE
     DEALLOCATE ( GRD%RO)
     DEALLOCATE ( GRD%RO_AD)
     DEALLOCATE ( ROS%EVCR, ROS%EVAR )
     DEALLOCATE ( ROS%EVC, ROS%EVA )
! CONTROL STRUCTURE
     DEALLOCATE( CTL%NBD, CTL%IWA)
     DEALLOCATE( CTL%X_C, CTL%G_C)
     DEALLOCATE( CTL%L_C, CTL%U_C)
     DEALLOCATE( CTL%WA, CTL%SG, CTL%SGO, CTL%YG, CTL%YGO)
     DEALLOCATE( CTL%WA3 )
     DEALLOCATE( CTL%WS, CTL%WY)
     DEALLOCATE( CTL%SY, CTL%SS, CTL%YY)
     DEALLOCATE( CTL%WT, CTL%WN, CTL%SND)
     DEALLOCATE( CTL%Z_C, CTL%R_C, CTL%D_C, CTL%T_C)
! BAROTROPIC MODEL
   IF(DRV%BMD(DRV%KTR) .EQ. 1) THEN
     DEALLOCATE ( BMD%ITR)
     DEALLOCATE ( BMD%MST, BMD%MSU, BMD%MSV)
     DEALLOCATE ( BMD%HGT, BMD%HGU, BMD%HGV)
     DEALLOCATE ( BMD%DXU, BMD%DXV)
     DEALLOCATE ( BMD%DYU, BMD%DYV)
     DEALLOCATE ( BMD%A1, BMD%A2, BMD%A3)
     DEALLOCATE ( BMD%A4, BMD%A0, BMD%A00)
     DEALLOCATE ( BMD%BX, BMD%BY)
     DEALLOCATE ( BMD%B_X, BMD%B_Y)
     DEALLOCATE ( BMD%DNS)
     DEALLOCATE ( BMD%BXBY, BMD%RGH)
     DEALLOCATE ( BMD%ETB, BMD%UB, BMD%VB)
     DEALLOCATE ( BMD%ETN, BMD%UN, BMD%VN)
     DEALLOCATE ( BMD%ETA, BMD%UA, BMD%VA)
     DEALLOCATE ( BMD%ETM, BMD%UM, BMD%VM)
     DEALLOCATE ( BMD%DIV, BMD%CU, BMD%CV)
     DEALLOCATE ( BMD%DUX, BMD%DUY)
     DEALLOCATE ( BMD%DVX, BMD%DVY)
     DEALLOCATE ( BMD%ETX, BMD%ETY)
   ENDIF


END SUBROUTINE SAV_ITR
