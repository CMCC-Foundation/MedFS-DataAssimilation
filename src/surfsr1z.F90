SUBROUTINE SURFSR1Z

!.. SETUP OF RECURSIVE FILTER ..!
!
!   THIS IMPLEMENTS THE 2010 SUMMER REVISION
!   OF THE RECURSIVE FILTER IMPLEMENTATION
!   WITH RECURSIVE FILTER COEFFICIENTS
!   FROM DISHOMOGENEOUS AND NON-ISOTROPIC
!   CORRELATION LENGTH-SCALES.
!
!   AUTHOR : ANDREA STORTO - 08-2010

USE SET_KND
USE RECFILTER
USE GRD_STR
USE OBS_STR, ONLY : SLA
USE MYFRTPROF, ONLY : MYFRTPROF_WALL
USE IOUNITS, ONLY : IOUNOUT,IOUNLOG
USE RUN, ONLY : LLVARMDT

IMPLICIT NONE
CALL MYFRTPROF_WALL('SURFSR1Z: RECURSIVE FILTER SETUP FOR SSH',0)

IF(SLA%NO .EQ. 0) THEN
   WRITE(IOUNOUT,*) ' *** NO SLA OBS, LLVARMDT IS IMPOSED FALSE'
   WRITE(IOUNLOG,*) ' *** NO SLA OBS, LLVARMDT IS IMPOSED FALSE'
   LLVARMDT = .FALSE.
   CALL MYFRTPROF_WALL('SURFSR1Z: RECURSIVE FILTER SETUP FOR SSH',1)
   RETURN
ENDIF

ALLOCATE(SR_SCXZ(GRD%IM,GRD%JM),&
       & SR_SCYZ(GRD%IM,GRD%JM),&
       & SR_GALXZ(GRD%IM,GRD%JM),&
       & SR_GALYZ(GRD%IM,GRD%JM),&
       & SR_GBTXZ(GRD%IM,GRD%JM),&
       & SR_GBTYZ(GRD%IM,GRD%JM) )

  !... READ COEFFS
  WRITE(IOUNOUT,'(2X,A)',ADVANCE='NO') 'READ_RFC: READ RF COEFFS...'
  CALL READ_RFCZ(GRD%IM,GRD%JM,SR_SCXZ,SR_SCYZ,SR_GALXZ,SR_GALYZ,&
  & SR_GBTXZ,SR_GBTYZ )
  WRITE(IOUNOUT,'(2X,A)') 'DONE'
  CALL FLUSH(IOUNLOG)
  CALL FLUSH(IOUNOUT)

  LLRFINIT = .TRUE.

  !... DEFINE EXTENDED AREA
  WRITE(IOUNOUT,'(2X,A)',ADVANCE='NO') 'SURF: EXTENDED GRID...'
  CALL DEFEXTGRDSRZ
  WRITE(IOUNOUT,'(2X,A)') 'DONE'
  CALL FLUSH(IOUNLOG)

  !... REARRANGE ALPHAS ACCORDING TO EXTENDED AREA
  WRITE(IOUNOUT,'(2X,A)',ADVANCE='NO') 'SURF: REARRANGE ALPHAS...'
  CALL REARRALPHASRZ
  WRITE(IOUNOUT,'(2X,A)') 'DONE'
  CALL FLUSH(IOUNLOG)

DEALLOCATE(SR_SCXZ,&
       & SR_SCYZ,&
       & SR_GALXZ,&
       & SR_GALYZ,&
       & SR_GBTXZ,&
       & SR_GBTYZ )

RETURN
CALL MYFRTPROF_WALL('SURFSR1Z: RECURSIVE FILTER SETUP FOR SSH',1)
END SUBROUTINE SURFSR1Z
