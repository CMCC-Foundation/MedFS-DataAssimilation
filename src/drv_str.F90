#ifndef USE_POINTERS
#define DYNMEM ALLOCATABLE
#else
#define DYNMEM POINTER
#endif

MODULE DRV_STR

!-----------------------------------------------------------------------
!                                                                      !
! STRUCTURE FOR THE DRIVER OF THE OUTER LOOP                           !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

 USE SET_KND

IMPLICIT NONE

INTEGER(I4), PARAMETER :: MGRIDS=6

PUBLIC

   TYPE DRV_T

        INTEGER(I4)           ::  NTR          ! NO. OF OUTER ITERATIONS
        INTEGER(I4)           ::  KTR          ! OUTER ITERATION
        INTEGER(I4)           ::  IM           ! DIMENSION OF THE COARSE GRID
        INTEGER(I4)           ::  JM           ! DIMENSION OF THE COARSE GRID
#ifndef NECSX
        INTEGER(I4), DYNMEM  ::  GRID(:)      ! GRID NUMBER FOR THE CURRENT ITERRATION
        CHARACTER(LEN=50), DYNMEM  ::  CGRID(:)     ! GRID NAME FOR THE CURRENT ITERRATION
        REAL(R8),    DYNMEM  ::  RATIO(:)     ! RATIO BETWEEN SUCCESSIVE GRIDS
        INTEGER(I4), DYNMEM  ::  MASK(:)      ! MASK USED FOR HORIZONTAL COVARIANCES
        INTEGER(I4), DYNMEM  ::  BMD(:)       ! 1 - RUN BAROTROPIC MODEL, ELSE - DO NOT RUN
        INTEGER(I4), DYNMEM  ::  DDA(:)       ! 1 - DIVERGENCE DAMPING IN ANALYSIS, ELSE NO FILTER
        INTEGER(I4), DYNMEM  ::  DDI(:)       ! 1 - DIVERGENCE DAMPING IN INITIALISATION, ELSE NO FILTER
        INTEGER(I4), DYNMEM  ::  BAL(:)       ! 1 - DIVERGENCE DAMPING IN INITIALISATION, ELSE NO FILTER
        INTEGER(I4), DYNMEM  ::  BA2(:)       ! 1 - DIVERGENCE DAMPING IN INITIALISATION, ELSE NO FILTER
#else
        INTEGER(I4) ::  GRID(MGRIDS)  ! GRID NUMBER FOR THE CURRENT ITERRATION
        CHARACTER(LEN=50)  ::  CGRID(MGRIDS)     ! GRID NAME FOR THE CURRENT ITERRATION
        REAL(R8)  ::  RATIO(MGRIDS)   ! RATIO BETWEEN SUCCESSIVE GRIDS
        INTEGER(I4) ::  MASK(MGRIDS)  ! MASK USED FOR HORIZONTAL COVARIANCES
        INTEGER(I4) ::  BMD(MGRIDS)   ! 1 - RUN BAROTROPIC MODEL, ELSE - DO NOT RUN
        INTEGER(I4) ::  DDA(MGRIDS)   ! 1 - DIVERGENCE DAMPING IN ANALYSIS, ELSE NO FILTER
        INTEGER(I4) ::  DDI(MGRIDS)   ! 1 - DIVERGENCE DAMPING IN INITIALISATION, ELSE NO FILTER
        INTEGER(I4) ::  BAL(MGRIDS)   ! 1 - DIVERGENCE DAMPING IN INITIALISATION, ELSE NO FILTER
        INTEGER(I4) ::  BA2(MGRIDS)   ! 1 - DIVERGENCE DAMPING IN INITIALISATION, ELSE NO FILTER
#endif
        REAL(R8)              ::  F_CI         ! INITAL COST FUNCTION
        REAL(R8),    DYNMEM  ::  RO   (:,:,:) ! VECTOR V
        REAL(R8)              ::  F_C          ! COST FUNCTION
        REAL(R8),    DYNMEM  ::  RO_AD(:,:,:) ! OBSERVATIONAL PART OF THE COST FUNCTION GRADIENT
        REAL(R8),    DYNMEM  ::    MSK(:,:)   ! MASK OF THE OLD GRID

   END TYPE DRV_T

   TYPE (DRV_T)              :: DRV

END MODULE DRV_STR
