#ifndef USE_POINTERS
#define DYNMEM ALLOCATABLE
#else
#define DYNMEM POINTER
#endif

MODULE CTL_STR

!-----------------------------------------------------------------------
!                                                                      !
! COST FUNCTION, CONTROL VECTOR AND OPTIMISATION ARRAYS                !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

 USE SET_KND

IMPLICIT NONE

PUBLIC

! ---
! STRUCTURE FOR LBFGS

   TYPE LBFGS_T

        CHARACTER*10              ::  ALG        ! MINIMZER (M1QN3 VS L-BFGS-B)
        INTEGER(I4)               ::  NIMAX      ! MAX NO ITERATIONS  (M1QN3)
        INTEGER(I4)               ::  NSMAX      ! MAX NO SIMULATIONS (M1QN3)
        INTEGER(I4)               ::  N          ! SIZE OF THE OPTIMISATION VECTOR
        INTEGER(I4)               ::  M          ! NUMBER OF COPIES TO BE SAVED
        INTEGER(I4)               ::  KITER      ! ITERATION INDEX
        CHARACTER*60              ::  TASK
        CHARACTER*60              ::  CSAVE
        LOGICAL, DIMENSION(4)     ::  LSAVE
        INTEGER(I4), DIMENSION(44)::  ISAVE
        INTEGER(I4), DYNMEM      ::  NBD(:), IWA(:)
        INTEGER(I4)               ::  IPRINT
        INTEGER(I4)               ::  IREST
        INTEGER(I4)               ::  IFLAG
        INTEGER(I4)               ::  GIPRINT(2)
        INTEGER(I4)               ::  GMETHOD
        REAL(R8)                  ::  F_B        ! THE BACKGROUND COST FUNCTION
        REAL(R8)                  ::  F_O        ! THE OBSERVATIONAL COST FUNCTION
        REAL(R8)                  ::  F_M        ! THE MDT COST FUNCTION
        REAL(R8)                  ::  F_E        ! THE MODEL ERRROR COST FUNCTION
        REAL(R8)                  ::  F_H        ! THE HYBRID ERROR COST FUNCTION
        REAL(R8)                  ::  F_A        ! THE ALPHA WEIGTH COST FUNCTION
        REAL(R8)                  ::  F_3        ! THE 3D EOFS COST FUNCTION
        REAL(R8)                  ::  F_T        ! THE TOA COST FUNCTION
#ifdef DOUBLE_INSTEAD_OF_R8
        DOUBLE PRECISION          ::  F_C, FACTR ! THE COST FUNCTION, ACCURACY
        DOUBLE PRECISION          ::  PGTOL, PGPER ! STOPPING CRITERIA, PERCENTAGE OF INITIAL GRADIENT
        DOUBLE PRECISION,  &
                 DIMENSION(29)    ::  DSAVE
        DOUBLE PRECISION,  &
                 DYNMEM          ::  X_C(:)     ! THE CONTROL VECTOR (BACKGROUND - ANALYSES)
        DOUBLE PRECISION,  &
                 DYNMEM          ::  G_C(:)     ! THE GRADIENT OF F_C
        DOUBLE PRECISION,  &
                 DYNMEM          ::  L_C(:), U_C(:)
        DOUBLE PRECISION,  &
                 DYNMEM          ::  WA(:), SG(:), SGO(:), YG(:), YGO(:),WA3(:),&
                                  WS(:,:), WY(:,:), SY(:,:), SS(:,:),       &
                                  YY(:,:), WT(:,:), WN(:,:), SND(:,:),      &
                                  Z_C(:), R_C(:), D_C(:), T_C(:)        ! WORKING ARRAYS
#else
        REAL(R8)          ::  F_C, FACTR,EPS ! THE COST FUNCTION, ACCURACY
        REAL(R8)          ::  PGTOL, PGPER ! STOPPING CRITERIA, PERCENTAGE OF INITIAL GRADIENT
        REAL(R8),  &
                 DIMENSION(29)    ::  DSAVE
        REAL(R8),  &
                 DYNMEM          ::  X_C(:)     ! THE CONTROL VECTOR (BACKGROUND - ANALYSES)
        REAL(R8),  &
                 DYNMEM          ::  G_C(:)     ! THE GRADIENT OF F_C
        REAL(R8),  &
                 DYNMEM          ::  L_C(:), U_C(:)
        REAL(R8),  &
                 DYNMEM          ::  WA(:), SG(:), SGO(:), YG(:), YGO(:),WA3(:), &
                                  WS(:,:), WY(:,:), SY(:,:), SS(:,:),       &
                                  YY(:,:), WT(:,:), WN(:,:), SND(:,:),      &
                                  Z_C(:), R_C(:), D_C(:), T_C(:)        ! WORKING ARRAYS
        REAL(R8),  DYNMEM        ::  GOLD(:),GW(:)
        LOGICAL                  ::  FINISH
#endif


   END TYPE LBFGS_T

   TYPE (LBFGS_T)                 :: CTL

   LOGICAL           :: LL_BFEND


END MODULE CTL_STR
