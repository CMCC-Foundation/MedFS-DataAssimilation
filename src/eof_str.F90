#ifndef USE_POINTERS
#define DYNMEM ALLOCATABLE
#else
#define DYNMEM POINTER
#endif

MODULE EOF_STR

!-----------------------------------------------------------------------
!                                                                      !
! STRUCTURE OF EOFS                                                    !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
! VERSION 2: A.STORTO   2010                                           !
!-----------------------------------------------------------------------

 USE SET_KND

IMPLICIT NONE

PUBLIC

   INTEGER(I4)  ::  NPRINTEOF
   LOGICAL :: LSQUAREDEVA
   LOGICAL :: LL_MDEOFS
   LOGICAL :: LL_WRITE_MDEOFS = .FALSE.
   LOGICAL :: LLFLOWDEP       = .FALSE.
   LOGICAL :: LLFLOWDEP_CALIBR= .FALSE.

   TYPE EOF_T

        CHARACTER(LEN=300)    ::  EOF_FILE     ! FILE CONTAINING EOFS
        INTEGER(I4)           ::  NEOF         ! NO. OF EOFS
        INTEGER(I4)           ::  NREG         ! NO. OF REGIONS
        INTEGER(I4)           ::  KMT          ! NO. OF LEVELS OF EOFS
        REAL(R8) :: DEP_MAX
        REAL(R8),    DYNMEM  ::  EVCR(:,:,:)  ! EIGENVECTORS ON REGIONS
        REAL(R8),    DYNMEM  ::  EVAR(:,:)    ! EIGENVALUES ON REGIONS
#ifdef opt_huge_memory
        REAL(R8),    DYNMEM  ::  EVC(:,:,:,:) ! EIGENVECTORS
        REAL(R8),    DYNMEM  ::  EVA(:,:,:)   ! EIGENVALUES
#else
        REAL(R8),    DYNMEM  ::  EVC(:,:,:)   ! EIGENVECTORS
        REAL(R8),    DYNMEM  ::  EVA(:,:)     ! EIGENVALUES
#endif
        INTEGER(I4)    ::  EOGNX,EOGNY,EOGNZ  ! DIMENSIONS OF EOF
                                              ! ORIGINALGRID

   END TYPE EOF_T

   TYPE (EOF_T)                 :: ROS
   TYPE (EOF_T)                 :: ROS2

END MODULE EOF_STR
