SUBROUTINE VER_HOR202

!-----------------------------------------------------------------------
!                                                                      !
! APPLY HORIZONTAL FILTER                                              !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
! VERSION 2: S.DOBRICIC 2007                                           !
!     SYMMETRIC CALCULATION IN PRESENCE OF COASTAL BOUNDARIES          !
!     ETA_AD, TEM_AD, AND SAL_AD ARE HERE TEMPORARY ARRAYS             !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE EOF_STR
 USE RECFILTER, ONLY : SCXNS,SCYNS,FCT,LLTSAPART,SCXNZ,SCYNZ
 USE DRV_STR
 USE CTL_STR
 USE IOUNITS
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL
 USE RUN, ONLY : LLVARMDT
 USE RECFILTER, ONLY : LLRFSR
 USE DEBUGTMP
 USE LBCLNK

 IMPLICIT NONE

 REAL(R8), DIMENSION (GRD%IM,GRD%JM)  :: UD, VD
 INTEGER(I4)    :: I,J,K,KL(2*GRD%KM),K1
 CHARACTER(LEN=20) :: CFILE

CALL MYFRTPROF_WALL('VER_HOR202: CONTROL TO PHYSICAL SPACE',0)

! ---
! VERTICAL EOFS
    GRD%CMDT(:,:) = GRD%CMDT(:,:) * GRD%CMDT_STDEV(:,:)

! ---
! LOAD TEMPORARY ARRAYS
    IF(DRV%MASK(DRV%KTR) .GT. 1) GRD%CMDT_AD(:,:  )   = GRD%CMDT(:,:  )

! ---
! X DIRECTION
    CALL RCFL_XZ( GRD%IM,GRD%JM,GRD%CMDT)

    IF( LL_LBCLNK ) CALL LBC_LNK( GRD%CMDT )

! ---
! SCALE BY THE SCALING FACTOR
       GRD%CMDT(:,:)   = GRD%CMDT(:,:)   * SCXNZ(:,:)

! ---
! Y DIRECTION
       CALL RCFL_YZ(GRD%IM,GRD%JM,GRD%CMDT)

! ---
! SCALE BY THE SCALING FACTOR
       GRD%CMDT(:,:)   = GRD%CMDT(:,:)   * SCYNZ(:,:)

! ---
! TRANSPOSE CALCULATION IN THE PRESENSE OF COASTAL BOUNDARIES
 IF(DRV%MASK(DRV%KTR) .GT. 1) THEN

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
       CALL RCFL_XZ_AD(GRD%IM,GRD%JM,GRD%CMDT_AD)

    IF( LL_LBCLNK ) CALL LBC_LNK( GRD%CMDT_AD )

! ---
! AVERAGE
        GRD%CMDT(:,:  )   = (GRD%CMDT(:,:  ) + GRD%CMDT_AD(:,:  ) ) * 0.5
ENDIF

! ---
! SCALE FOR BOUNDARIES
       GRD%CMDT(:,:)   = GRD%CMDT(:,:)   * FCT(:,:,1)

CALL FLUSH(IOUNLOG)

CALL MYFRTPROF_WALL('VER_HOR202: CONTROL TO PHYSICAL SPACE',1)
END SUBROUTINE VER_HOR202
