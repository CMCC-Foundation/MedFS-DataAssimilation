SUBROUTINE VER_HOR_AD202

!-----------------------------------------------------------------------
!                                                                      !
! APPLY VERTICAL EOFS AND HORIZONTAL FILTER - ADJOINT                  !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
! VERSION 2: S.DOBRICIC 2007                                           !
!     SYMMETRIC CALCULATION IN PRESENCE OF COASTAL BOUNDARIES          !
!     ETA, TEM, AND SAL ARE HERE TEMPORARY ARRAYS                      !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE EOF_STR
 USE DRV_STR
 USE CTL_STR
 USE RECFILTER, ONLY : SCXNS,SCYNS,FCT,LLTSAPART,LLRFSR,SCXNZ,SCYNZ
 USE RUN, ONLY : LLVARMDT
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL
 USE DEBUGTMP
 USE LBCLNK

 IMPLICIT NONE

 INTEGER(I4)    :: K,KL(2*GRD%KM)
 CHARACTER(LEN=20) :: CFILE

CALL MYFRTPROF_WALL('VER_HOR_AD202: ADJ OF CONTROL TO PHYSICAL SPACE',0)

! ---
! SCALE FOR BOUNDARIES
       GRD%CMDT_AD(:,:)   = GRD%CMDT_AD(:,:)   * FCT(:,:,1)

 IF(DRV%MASK(DRV%KTR) .GT. 1) THEN

! ---
! LOAD TEMPORARY ARRAYS
          GRD%CMDT(:,:  )    = GRD%CMDT_AD(:,:  )

! ---
! X DIRECTION
       CALL RCFL_XZ( GRD%IM, GRD%JM, GRD%CMDT)

    IF( LL_LBCLNK ) CALL LBC_LNK( GRD%CMDT )

! ---
! SCALE BY THE SCALING FACTOR
       GRD%CMDT(:,:)   = GRD%CMDT(:,:)   * SCXNZ(:,:)

! ---
! Y DIRECTION
       CALL RCFL_YZ(GRD%IM, GRD%JM,GRD%CMDT)

! ---
! SCALE BY THE SCALING FACTOR
       GRD%CMDT(:,:)   = GRD%CMDT(:,:)   * SCYNZ(:,:)

 ENDIF

! ---
! SCALE BY THE SCALING FACTOR
       GRD%CMDT_AD(:,:)   = GRD%CMDT_AD(:,:)   * SCYNZ(:,:)

! ---
! Y DIRECTION
       CALL RCFL_YZ_AD( GRD%IM, GRD%JM,GRD%CMDT_AD)

! ---
! SCALE BY THE SCALING FACTOR
       GRD%CMDT_AD(:,:)   = GRD%CMDT_AD(:,:)   * SCXNZ(:,:)

! ---
! X DIRECTION
       CALL RCFL_XZ_AD( GRD%IM, GRD%JM,GRD%CMDT_AD)
    IF( LL_LBCLNK ) CALL LBC_LNK( GRD%CMDT_AD )

! ---
! AVERAGE
 IF(DRV%MASK(DRV%KTR) .GT. 1) THEN
        GRD%CMDT_AD(:,:  )   = (GRD%CMDT_AD(:,:  ) + GRD%CMDT(:,:  ) ) * 0.5
 ENDIF

! ---
! VERTICAL EOFS
    GRD%CMDT_AD(:,:) = GRD%CMDT_AD(:,:) * GRD%CMDT_STDEV(:,:)

 CALL MYFRTPROF_WALL('VER_HOR_AD202: ADJ OF CONTROL TO PHYSICAL SPACE',1)

END SUBROUTINE VER_HOR_AD202
