#ifndef USE_POINTERS
#define DYNMEM ALLOCATABLE
#else
#define DYNMEM POINTER
#endif

MODULE GRD_STR

!-----------------------------------------------------------------------
!                                                                      !
! STRUCTURE OF THE GRID                                                !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!
! A.S. : NOW CONTROL VECTORS AND THEIR ADJOINTS IN PHYSICAL SPACE
! ARE JUST POINTERS TO A GLOBAL ARRAY, IN ORDER TO OPTIMIZE THE
! OPENMP PARALLELIZATION AND PERMIT INDEPENDENT TREATMENT OF T AND S
! IN THE RECURSIVE FILTER ROUTINES.
!-----------------------------------------------------------------------

 USE SET_KND

IMPLICIT NONE

PUBLIC

   REAL(R8), ALLOCATABLE, TARGET  :: PSV(:,:,:) ! PHYSICAL SPACE VECTOR
   REAL(R8), ALLOCATABLE, TARGET  :: PSV_AD(:,:,:) ! ADJ OF PHYS SP VECTOR
   LOGICAL :: LL_TS, LL_UV, LL_SSH, LL_SST, LL_SSS, LL_TQ2, LL_ICE
   LOGICAL :: LL_TQ2_ASSIM, LL_ICE_ASSIM
   INTEGER(I4) :: PSV3D, PSV2D
   REAL(R8), ALLOCATABLE  :: PSVF(:,:,:) ! PHYSICAL SPACE VECTOR
   REAL(R8), ALLOCATABLE  :: PSVF_AD(:,:,:) ! ADJ OF PHYS SP VECTOR

   REAL(R8) :: MIMA_LL(4)

   TYPE GRID_T

        INTEGER(I4)              ::  I0           ! START X
        INTEGER(I4)              ::  J0           ! START Y
        INTEGER(I4)              ::  IM           ! NO. POINTS IN X DIRECTION
        INTEGER(I4)              ::  JM           ! NO. POINTS IN Y DIRECTION
        INTEGER(I4)              ::  KM           ! NO. POINTS IN Z DIRECTION
        INTEGER(I4)              ::  NPS          ! NO. OF OCEAN POINTS

        REAL(R8)                 ::  DLN          ! RESOLUTION IN THE X DIRECTION
        REAL(R8)                 ::  DLT          ! RESOLUTION IN THE Y DIRECTION

        REAL(R8),    DYNMEM     ::  RO(:,:,:)     ! REDUCED ORDER CONTROL VECTOR
        INTEGER(I4), DYNMEM     ::  REG(:,:)      ! MASK FOR EOF REGIONS
        REAL(R8),    DYNMEM     ::  MSK(:,:,:)    ! SEA-LAND MASK
        REAL(R8),    DYNMEM     ::  MVLOC(:,:,:)  ! MLD VERT LOC
        REAL(R8),    DYNMEM     ::  HGT(:,:)      ! TOPOGRAPHY
        REAL(R8),    DYNMEM     ::  ICE(:,:)      ! SEA-ICE FRACTION 0-1
        REAL(R8),    DYNMEM     ::  DISTC(:,:)    ! DISTANCE FROM COAST
        REAL(R8),    DYNMEM     ::    F(:,:)      ! CORIOLIS TERM

        REAL(R8),    DYNMEM     :: CMDT(:,:)      ! MDT INCREMENT AND ADJOINT
        REAL(R8),    DYNMEM     :: CMDT_AD(:,:)   ! FOR VARMDT
        REAL(R8),    DYNMEM     :: CMDT_STDEV(:,:)! MDT ERROR

!... USE ONE HUGE CONTROL VECTOR IN PHYSICAL SPACE
!    AND MAKE THE INDIVIDUAL VARS POINT TO IT
!    POSSIBILITY TO USE ALLOCATABLE INDIVIDUAL MATRICES
!    IS RETAINED THROUGH PP KEY
        REAL(R8),    POINTER     ::  TEM(:,:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  SAL(:,:,:)   ! SALINITY INCREMENT
        REAL(R8),    POINTER     ::  UVL(:,:,:)   ! U COMP
        REAL(R8),    POINTER     ::  VVL(:,:,:)   ! V COMP
        REAL(R8),    POINTER     ::  ETA(:,:)     ! SEA LEVEL INCREMENT
        REAL(R8),    POINTER     ::  SST(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  SSS(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  T2M(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  Q2M(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  SIC(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  SIT(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  ICU(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  ICV(:,:)    ! TEMPERATURE INCREMENT

        REAL(R8),    POINTER     ::  TEM_AD(:,:,:)   ! TEMPERATURE ADJOINT
        REAL(R8),    POINTER     ::  SAL_AD(:,:,:)   ! SALINITY ADJOINT
        REAL(R8),    POINTER     ::  UVL_AD(:,:,:)   ! U COMP ADJOINT
        REAL(R8),    POINTER     ::  VVL_AD(:,:,:)   ! V COMP ADJOINT
        REAL(R8),    POINTER     ::  ETA_AD(:,:)     ! SEA LEVEL ADJOINT
        REAL(R8),    POINTER     ::  SST_AD(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  SSS_AD(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  T2M_AD(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  Q2M_AD(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  SIC_AD(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  SIT_AD(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  ICU_AD(:,:)    ! TEMPERATURE INCREMENT
        REAL(R8),    POINTER     ::  ICV_AD(:,:)    ! TEMPERATURE INCREMENT

        REAL(R8),    DYNMEM     ::  TEMB(:,:,:,:) ! TEMPERATURE BACKGROUND
        REAL(R8),    DYNMEM     ::  SALB(:,:,:,:) ! SALINITY BACKGROUND
        REAL(R8),    DYNMEM     ::  UVLB(:,:,:)   ! U COMPONNET OF VELOCITY BACKGROUND
        REAL(R8),    DYNMEM     ::  VVLB(:,:,:)   ! V COMPONNET OF VELOCITY BACKGROUND
        REAL(R8),    DYNMEM     ::  ETAB(:,:,:)   ! SEA LEVEL BACKGROUND
        REAL(R8),    DYNMEM     ::  T2MB(:,:,:)   ! SEA LEVEL BACKGROUND
        REAL(R8),    DYNMEM     ::  Q2MB(:,:,:)   ! SEA LEVEL BACKGROUND
        REAL(R8),    DYNMEM     ::  SICB(:,:,:)   ! SEA LEVEL BACKGROUND
        REAL(R8),    DYNMEM     ::  SITB(:,:,:)   ! SEA LEVEL BACKGROUND
        REAL(R8),    DYNMEM     ::  ICUB(:,:,:)   ! SEA LEVEL BACKGROUND
        REAL(R8),    DYNMEM     ::  ICVB(:,:,:)   ! SEA LEVEL BACKGROUND

        REAL(R8),    DYNMEM     ::  MDT(:,:)      ! MEAN DYNAMIC TOPOGRAPHY
        REAL(R8),    DYNMEM     ::  SLA(:,:)      ! SEA LEVEL ANOMALY

        REAL(R8),    DYNMEM     ::  RO_AD(:,:,:)    ! REDUCED ORDER CONTROL VECTOR ADJOINT
        REAL(R8),    DYNMEM     ::  DNS(:,:,:)      ! DENSITY
        REAL(R8),    DYNMEM     ::  B_X(:,:,:)      ! BOUYANCY FORCE
        REAL(R8),    DYNMEM     ::  B_Y(:,:,:)      ! BOUYANCY FORCE
        REAL(R8),    DYNMEM     ::  BX(:,:)         ! BOUYANCY FORCE INTEGRAL
        REAL(R8),    DYNMEM     ::  BY(:,:)         ! BOUYANCY FORCE INTEGRAL

        REAL(R8),    DYNMEM     ::  LON(:,:)       ! LONGITUDE
        REAL(R8),    DYNMEM     ::  LAT(:,:)       ! LATITUDE
        REAL(R8),    DYNMEM     ::  DEP(:)       ! DEPTH

        REAL(R8),    DYNMEM     ::  HVST(:,:,:)  ! HORIZONTAL DIFF. COEF. FOR TEMPERATURE
        REAL(R8),    DYNMEM     ::  HVSS(:,:,:)  ! HORIZONTAL DIFF. COEF. FOR SALINITY
        REAL(R8),    DYNMEM     ::  HVSP(:,:,:)  ! HORIZONTAL DIFF. COEF. FOR STREAM FUNCTION
        REAL(R8),    DYNMEM     ::  VVST(:,:,:)  ! VERTICAL DIFF. COEF. FOR TEMPERATURE
        REAL(R8),    DYNMEM     ::  VVSS(:,:,:)  ! VERTICAL DIFF. COEF. FOR SALINITY
        REAL(R8),    DYNMEM     ::  VVSP(:,:,:)  ! VERTICAL DIFF. COEF. FOR STREAM FUNCTION
        REAL(R8),    DYNMEM     ::  DX(:,:)      ! DX
        REAL(R8),    DYNMEM     ::  DY(:,:)      ! DY
        REAL(R8),    DYNMEM     ::  DZ(:)        ! DZ
        REAL(R8),    DYNMEM     ::  DXDY(:,:)    ! DX*DY
        REAL(R8)                ::  ADXDY       ! MEAN DX*DY

        REAL(R8),    DYNMEM     ::  MSR(:,:,:)   ! SEA-LAND MASK USED IN THE RECURSIVE FILTER

        ! HYBRID FORMULATION
        REAL(R8),    DYNMEM     ::  ROH(:,:,:)     ! REDUCED ORDER CONTROL VECTOR
        REAL(R8),    DYNMEM     ::  ROH_AD(:,:,:)    ! REDUCED ORDER CONTROL VECTOR ADJOINT
        REAL(R8),    DYNMEM     ::  MVLOC_H(:,:,:)  ! MLD VERT LOC

   END TYPE GRID_T

   TYPE (GRID_T)                 :: GRD

END MODULE GRD_STR
