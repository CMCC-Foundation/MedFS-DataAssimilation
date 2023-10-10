#ifndef USE_POINTERS
#define DYNMEM ALLOCATABLE
#else
#define DYNMEM POINTER
#endif

MODULE BMD_STR

!-----------------------------------------------------------------------
!                                                                      !
! STRUCTURE FOR THE BAROTROPIC MODEL                                   !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------

 USE SET_KND

IMPLICIT NONE

PUBLIC

   TYPE BMD_T

        INTEGER(I4)              ::  NCNT         ! MAXIMUM NUMBER OF ITERATIONS IN THE IMPLICIT SOLVER
        REAL(R8)                 ::  OVR          ! OVER-RELAXATION FACTOR
        REAL(R8)                 ::  RESEM        ! STOPPING CRITERIA
        REAL(R8)                 ::  BNM          ! NUMBER OF SEA POINTS

        REAL(R8)                 ::  G            ! GRAVIATIONAL ACCELERATION
        REAL(R8)                 ::  DT           ! TIME STEP
        INTEGER(I4)              ::  NSTP         ! NUMBER OF TIME STEPS PER DAY
        REAL(R8)                 ::  NDY          ! NUMBER OF SIMULATION DAYS
        REAL(R8)                 ::  ADY          ! NUMBER OF AVERAGING DAYS
        INTEGER(I4)              ::  NSTPS        ! NUMBER OF TIME STEPS OF THE MAIN LOOP
        INTEGER(I4)              ::  NSTPA        ! NUMBER OF TIME STEPS FOR AVERAGING
        REAL(R8)                 ::  ALP1         ! WEIGHTING FACTOR IN THE TRAPEZOIDAL SCHEME
        REAL(R8)                 ::  ALP2         ! WEIGHTING FACTOR IN THE TRAPEZOIDAL SCHEME
        REAL(R8)                 ::  FC1          ! FRICTION INTENSITY
        REAL(R8)                 ::  FC2          ! FRICTION INTENSITY
        REAL(R8)                 ::  DF1          ! FRICTION INTENSITY
        REAL(R8)                 ::  DF2          ! FRICTION INTENSITY

        INTEGER(I4), DYNMEM     ::  ITR(:)       ! NUMBER OF ITERATIONS IN THE SOLVER
        REAL(R8),    DYNMEM     ::  MST(:,:)     ! SEA-LAND MASK ON T POINTS
        REAL(R8),    DYNMEM     ::  MSU(:,:)     ! SEA-LAND MASK ON U POINTS
        REAL(R8),    DYNMEM     ::  MSV(:,:)     ! SEA-LAND MASK ON V POINTS
        REAL(R8),    DYNMEM     ::  HGT(:,:)     ! DEPTH ON T POINTS
        REAL(R8),    DYNMEM     ::  HGU(:,:)     ! DEPTH ON U POINTS
        REAL(R8),    DYNMEM     ::  HGV(:,:)     ! DEPTH ON V POINTS
        REAL(R8),    DYNMEM     ::  DXU(:,:)     ! DX ON U POINTS
        REAL(R8),    DYNMEM     ::  DYU(:,:)     ! DY ON U POINTS
        REAL(R8),    DYNMEM     ::  DXV(:,:)     ! DX ON V POINTS
        REAL(R8),    DYNMEM     ::  DYV(:,:)     ! DY ON V POINTS
        REAL(R8),    DYNMEM     ::   A1(:,:)     ! CONSTANT
        REAL(R8),    DYNMEM     ::   A2(:,:)     ! CONSTANT
        REAL(R8),    DYNMEM     ::   A3(:,:)     ! CONSTANT
        REAL(R8),    DYNMEM     ::   A4(:,:)     ! CONSTANT
        REAL(R8),    DYNMEM     ::   A0(:,:)     ! CONSTANT
        REAL(R8),    DYNMEM     ::  A00(:,:)     ! CONSTANT
        REAL(R8),    DYNMEM     ::   BX(:,:)     ! BOUYANCY GRADIENT IN X DIRECTION (VERT. INT.)
        REAL(R8),    DYNMEM     ::   BY(:,:)     ! BOUYANCY GRADIENT IN Y DIRECTION (VERT. INT.)
        REAL(R8),    DYNMEM     ::  B_X(:,:,:)   ! BOUYANCY GRADIENT IN X DIRECTION
        REAL(R8),    DYNMEM     ::  B_Y(:,:,:)   ! BOUYANCY GRADIENT IN Y DIRECTION
        REAL(R8),    DYNMEM     ::  DNS(:,:,:)   ! DENSITY
        REAL(R8),    DYNMEM     ::  BXBY(:,:)    !
        REAL(R8),    DYNMEM     ::   RGH(:,:)    !
        REAL(R8),    DYNMEM     ::   ETB(:,:)    ! ETA AT T-1
        REAL(R8),    DYNMEM     ::    UB(:,:)    ! U AT T-1
        REAL(R8),    DYNMEM     ::    VB(:,:)    ! V AT T-1
        REAL(R8),    DYNMEM     ::   ETN(:,:)    ! ETA AT T
        REAL(R8),    DYNMEM     ::    UN(:,:)    ! U AT T
        REAL(R8),    DYNMEM     ::    VN(:,:)    ! V AT T
        REAL(R8),    DYNMEM     ::   ETA(:,:)    ! ETA AT T+1
        REAL(R8),    DYNMEM     ::    UA(:,:)    ! U AT T+1
        REAL(R8),    DYNMEM     ::    VA(:,:)    ! V AT T+1
        REAL(R8),    DYNMEM     ::   ETM(:,:)    ! AVERAGED ETA
        REAL(R8),    DYNMEM     ::    UM(:,:)    ! AVERAGED U
        REAL(R8),    DYNMEM     ::    VM(:,:)    ! AVERAGED V
        REAL(R8),    DYNMEM     ::   DIV(:,:)    ! DIVERGENCE AT T-1
        REAL(R8),    DYNMEM     ::    CU(:,:)    ! CORIOLIS TERM ON U POINTS
        REAL(R8),    DYNMEM     ::    CV(:,:)    ! CORIOLIS TERM ON V POINTS
        REAL(R8),    DYNMEM     ::   DUX(:,:)    ! FRICTION ON U
        REAL(R8),    DYNMEM     ::   DUY(:,:)    ! FRICTION ON U
        REAL(R8),    DYNMEM     ::   DVX(:,:)    ! FRICTION ON V
        REAL(R8),    DYNMEM     ::   DVY(:,:)    ! FRICTION ON V
        REAL(R8),    DYNMEM     ::   ETX(:,:)    ! FREE SURFACE GRADIENT AT T-1
        REAL(R8),    DYNMEM     ::   ETY(:,:)    ! FREE SURFACE GRADIENT AT T-1
        REAL(R8),    DYNMEM     ::   BFU(:,:)    ! BOTTOM FRICTION
        REAL(R8),    DYNMEM     ::   BFV(:,:)    ! BOTTOM FRICTION


   END TYPE BMD_T

   TYPE (BMD_T)                 :: BMD

END MODULE BMD_STR
