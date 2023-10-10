MODULE MPIREL

!! MPI RELATED DEFINITIONS
!!
!!
!! THESE VARIABLES ARE SETUP IN <SUMPI>

USE SET_KND, ONLY : I4
#ifndef NOMPI
USE MPI
#endif

IMPLICIT NONE

PUBLIC
SAVE

!... MPI TO USE/NOT
LOGICAL :: LLUSEMPI
!... MPI USED/NOT
LOGICAL :: LMPION
!... MPI WE EA BOUNDARY CONNECTION
LOGICAL :: EA_BOUNDARY,WE_BOUNDARY,SO_BOUNDARY,NO_BOUNDARY
!... RANK OF PROCS
INTEGER(KIND=I4) :: NPROCS
!... WHICH PROCESSOR AM I? STARTING FROM 1 TO NPROCS
INTEGER(KIND=I4) :: MYPROC
!... PROCESSOR DEDICATED FOR I/O UNITS
INTEGER(KIND=I4) :: NIOMASTER
!... PROCESSOR DEDICATED FOR I/O UNITS, AUXILIARY
INTEGER(KIND=I4) :: NIOMASTERB
!... MPI DOMAIN
CHARACTER(LEN=4) :: CMPIDOM
!... MPI DOMAIN LIMITS
INTEGER(KIND=I4), ALLOCATABLE :: MP_SD(:), MP_SDST(:), MP_XST(:), MP_XEN(:), MP_YST(:), MP_YEN(:)
!... MPI DOMAIN MAP
INTEGER(KIND=I4), ALLOCATABLE :: PROC_XMED(:), PROC_YMED(:), PROC_XMED_ORDERED(:), PROC_YMED_ORDERED(:), PROC_WE(:), &
& PROC_EA(:),PROC_NO(:),PROC_SO(:),PROC_POSX_IN_MAP(:),PROC_POSY_IN_MAP(:),MAP_XMED(:), &
&  MAP_YMED(:), MPI_MAP(:,:)
INTEGER(KIND=I4), ALLOCATABLE :: LOCAL_OVERLAP_WE_XST(:), LOCAL_OVERLAP_WE_XEN(:), LOCAL_OVERLAP_EA_XST(:), &
& LOCAL_OVERLAP_EA_XEN(:),   LOCAL_OVERLAP_SO_YST(:), LOCAL_OVERLAP_SO_YEN(:),LOCAL_OVERLAP_NO_YST(:),LOCAL_OVERLAP_NO_YEN(:)
! NEXT ARE NECESSARY IN CASE OF GRID EXTENSION IN EA-WE
INTEGER(KIND=I4), ALLOCATABLE :: LOCAL_OVERLAP_SO_XST(:), LOCAL_OVERLAP_SO_XEN(:),LOCAL_OVERLAP_NO_XST(:),LOCAL_OVERLAP_NO_XEN(:)
LOGICAL, ALLOCATABLE, DIMENSION(:)  :: WE_BOUND_ARRAY,EA_BOUND_ARRAY,NO_BOUND_ARRAY,SO_BOUND_ARRAY
INTEGER(KIND=I4) :: NXDOM,NYDOM

#ifndef NOMPI
INTERFACE mppsum
     MODULE PROCEDURE mppsum_int, mppsum_real4, mppsum_real8
END INTERFACE
#else
INTERFACE mppsum
     MODULE PROCEDURE sum_int, sum_real4, sum_real8
END INTERFACE
#endif

CONTAINS

   SUBROUTINE sum_int(ktab)
    INTEGER(I4), INTENT(inout) ::   ktab
   END SUBROUTINE sum_int

   SUBROUTINE sum_real8(ptab)
    REAL(8), INTENT(inout)           ::   ptab   ! input scalar
   END SUBROUTINE sum_real8

   SUBROUTINE sum_real4(ptab)
    REAL(4), INTENT(inout)           ::   ptab   ! input scalar
   END SUBROUTINE sum_real4

#ifndef NOMPI
   SUBROUTINE mppsum_int(ktab)
    INTEGER(I4), INTENT(inout) ::   ktab
    INTEGER(I4) :: ierror, iwork
    CALL mpi_allreduce(ktab,iwork,1,mpi_integer,mpi_sum,MPI_COMM_WORLD,ierror)
    ktab = iwork
   END SUBROUTINE mppsum_int

   SUBROUTINE mppsum_real8(ptab)
    REAL(8), INTENT(inout)           ::   ptab   ! input scalar
    INTEGER(I4)  ::   ierror
    REAL(8) ::   zwork
    CALL mpi_allreduce(ptab,zwork,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierror)
    ptab = zwork
   END SUBROUTINE mppsum_real8

   SUBROUTINE mppsum_real4(ptab)
    REAL(4), INTENT(inout)           ::   ptab   ! input scalar
    INTEGER(I4)  ::   ierror
    REAL(4) ::   zwork
    CALL mpi_allreduce(ptab,zwork,1,mpi_real,mpi_sum,MPI_COMM_WORLD,ierror)
    ptab = zwork
   END SUBROUTINE mppsum_real4
#endif

END MODULE MPIREL
