MODULE PROF_CORREL

USE SET_KND
USE OBS_STR
USE GRD_STR
USE IOUNITS
USE MPIREL
USE MYFRTPROF
USE RUN
USE OCEANTOOLS
USE MYNETCDF

IMPLICIT NONE

INTEGER(I4)       :: NN_PROF_CORREL = 0
LOGICAL           :: LL_PROFCORR_ONLYP = .FALSE.
REAL(R8), ALLOCATABLE :: PMAT(:,:,:)
INTEGER(I4), PARAMETER :: MAXPROFS=9999
INTEGER(I4) :: NPROFS, NOBSPT(MAXPROFS), &
& JOSTART(MAXPROFS), JOEND(MAXPROFS)
REAL(R8) :: RDSC(MAXPROFS)

REAL(R8), PARAMETER :: DEP_MAX = 500._R8

CONTAINS

SUBROUTINE JOINS_CORREL(JOS)

IMPLICIT NONE

REAL(R8), INTENT(OUT) :: JOS
INTEGER(I4) :: JPRF, IMAX, OFF
REAL(R8), ALLOCATABLE :: ZTMP(:,:)
REAL(R8) :: ZT(1,1)

IMAX= MAXVAL( NOBSPT )
ALLOCATE(ZTMP(IMAX,1))

WRITE(1591+MYPROC,*)
OFF = SLA%NO

JOS=0._R8
DO JPRF=1,NPROFS

   ZTMP(1:NOBSPT(JPRF),1) = OBS%INC((OFF+JOSTART(JPRF)):(OFF+JOEND(JPRF)))-&
   & OBS%RES((OFF+JOSTART(JPRF)):(OFF+JOEND(JPRF)))

   ZT = MATMUL(TRANSPOSE(ZTMP(1:NOBSPT(JPRF),:)),&
   & MATMUL(PMAT(JPRF,1:NOBSPT(JPRF),1:NOBSPT(JPRF)),ZTMP(1:NOBSPT(JPRF),:)))

   JOS = JOS + 0.5_R8 * ZT(1,1)

ENDDO

WRITE(1591+MYPROC,*)
DEALLOCATE(ZTMP)

END SUBROUTINE JOINS_CORREL

SUBROUTINE GRAINS_CORREL(N,GRA)

IMPLICIT NONE

INTEGER(I4), INTENT(IN) :: N
REAL(R8), INTENT(OUT) :: GRA(N)
REAL(R8) :: GRAT(N,1)
INTEGER(I4) :: JPRF, IMAX, OFF
REAL(R8), ALLOCATABLE :: ZTMP(:,:)

IMAX= MAXVAL( NOBSPT )
ALLOCATE(ZTMP(IMAX,1))
OFF = SLA%NO

DO JPRF=1,NPROFS

   ZTMP(1:NOBSPT(JPRF),1) = OBS%INC((OFF+JOSTART(JPRF)):(OFF+JOEND(JPRF)))-&
   & OBS%RES((OFF+JOSTART(JPRF)):(OFF+JOEND(JPRF)))

   GRAT( JOSTART(JPRF):JOEND(JPRF),: ) = &
   & MATMUL(PMAT(JPRF,1:NOBSPT(JPRF),1:NOBSPT(JPRF)),ZTMP(1:NOBSPT(JPRF),:))

ENDDO

GRA(1:N)=GRAT(1:N,1)
DEALLOCATE(ZTMP)

END SUBROUTINE GRAINS_CORREL

SUBROUTINE PREP_INSRINV

!... SEARCH OUT INS ALONG-PROF CO-LOCATIONS
!    FOR ERROR CORRELATION COMPUTATION
!... A.S.

IMPLICIT NONE

REAL(R8) :: CR2, DIST, DET, ZCHK, ZCT, CRT, COR, MAD, MID
REAL(R8), ALLOCATABLE :: RTMP(:,:),RODENS(:),RMDENS(:,:,:),DPROF(:),PQ(:)
INTEGER(I4) :: JO, JPRF, KO, JO2, KO2
INTEGER(I4) :: MA, IERR, KPROF, KTYP
INTEGER(I4) :: JD, NDIST, IDN=0, JP, JK, NREJ
LOGICAL :: LLDEB
CHARACTER(LEN=20) :: CFMT

#include "obs_events.h"

 CALL MYFRTPROF_WALL('PREP_INSRINV: PREPARING FOR INS R INVERSION',0)

 WRITE(IOUNLOG,*)
 WRITE(IOUNLOG,*) ' *** PREPARING FOR INS R INVERSION'
 WRITE(IOUNLOG,*) 
 WRITE(IOUNLOG,*) ' NUMBER OF INS OBS ',INS%NO

 IF(LL_HUBERQC_INS) CALL ABOR1('INSRINV NOT SUPPORTED WITH HUBERQC, ABORTING')

 IF(NN_PROF_CORREL.NE. 1 ) THEN
   CALL ABOR1('INSRINV: NN_PROF_CORREL OPTION NOT SUPPORTED')
 ENDIF

 NPROFS = 0
 KPROF  = 0
 KTYP   = 0
 
 ! ALLOCATE / READ DENSITY
 ALLOCATE(RMDENS(GRD%IM,GRD%JM,GRD%KM))
 ALLOCATE(RODENS(INS%NO))

 ALLOCATE(PQ(NPQ),DPROF(GRD%KM))

 CALL READ_DENS(GRD%IM,GRD%JM,GRD%KM,RMDENS)
 IF( MYPROC.EQ.1 ) LLDEB =.TRUE.

 NREJ = 0

 LOOPINS : DO JO=1,INS%NO
    IF( INS%FLC(JO) .EQ. 0 ) CYCLE LOOPINS 
    ! DENSITY CALCUL.
    RODENS(JO) = 0._R8
    DO JP=1,NPQ
       RODENS(JO) = RODENS(JO) + &
          INS%PQ(JO,JP)    *RMDENS(INS%IB(JO,JP),INS%JB(JO,JP),INS%KB(JO)  ) + &
          INS%PQ(JO,JP+NPQ)*RMDENS(INS%IB(JO,JP),INS%JB(JO,JP),INS%KB(JO)+1)
    ENDDO
    IF( INS%OTYPE(JO) .EQ. KTYP .AND. INS%PROF(JO) .EQ. KPROF ) THEN
       NOBSPT( NPROFS ) = NOBSPT( NPROFS ) + 1
       CYCLE LOOPINS
    ENDIF
    DPROF  = -999._R8
    MAD    = 0._R8
    MID    = 2000._R8
    LOOPDEP : DO JK=1,GRD%KM-1
       IF( GRD%DEP(JK).GT.DEP_MAX) EXIT LOOPDEP
       DO JP=1,NPQ
          PQ(JP) = INS%PQ(JO,JP)*GRD%MSK(INS%IB(JO,JP),INS%JB(JO,JP),JK)
       ENDDO
       IF(SUM(PQ).GT.0._R8) THEN
         PQ = PQ / SUM(PQ)
         DPROF(JK) = 0._R8
         DO JP=1,NPQ
          DPROF(JK)=DPROF(JK)+PQ(JP)*RMDENS(INS%IB(JO,JP),INS%JB(JO,JP),JK)
         ENDDO
         IF(DPROF(JK).GT.MAD) MAD = DPROF(JK)
         IF(DPROF(JK).LT.MID) MID = DPROF(JK)
       ENDIF
    ENDDO LOOPDEP
    IF(MAD.EQ.0._R8.OR.MID.EQ.2000._R8) THEN
       WRITE(IOUNLOG,*) '/// I ', INS%IB(JO,:)
       WRITE(IOUNLOG,*) '/// J ', INS%JB(JO,:)
       WRITE(IOUNLOG,*) '/// P ', INS%PQ(JO,1:4)
       INS%FLC(JO) = 0
       INS%EVE(JO) = KEVE_INTE
       NREJ = NREJ + 1
       CYCLE LOOPINS
    ENDIF

    NPROFS = NPROFS + 1
    RDSC(NPROFS) = 0.25_R8 * (MAD-MID)
    IF(LLDEB) WRITE(IOUNLOG,*) JO,INS%IB(JO,1),INS%JB(JO,1),MID,MAD
    NOBSPT( NPROFS ) = 1
    KTYP=INS%OTYPE(JO)
    KPROF=INS%PROF(JO)
    JOSTART(NPROFS) = JO
    IF( NPROFS.GT.1) JOEND(NPROFS-1) = JO-1
 ENDDO LOOPINS
 JOEND(NPROFS)=INS%NO

 DEALLOCATE(RMDENS,PQ,DPROF)

 MA=MAXVAL( NOBSPT(1:NPROFS) )

 WRITE(IOUNLOG,*)
 WRITE(IOUNLOG,*) ' REPORT OF INS PROFS'
 WRITE(IOUNLOG,*) ' REJECTED  INS OBS   = ',NREJ
 WRITE(IOUNLOG,*) ' NUMBER OF INS PROFS = ',NPROFS
 DO JPRF=1,NPROFS
   WRITE(IOUNLOG,*) ' PROF NO',JPRF,'  -  ',NOBSPT(JPRF),JOSTART(JPRF),JOEND(JPRF),RDSC(JPRF)
 ENDDO

 WRITE(IOUNLOG,*) ' ALLOCATING PMAT :', NPROFS, MA, MA
 WRITE(IOUNLOG,*) ' ABOUT ', NINT(NPROFS*MA*MA*8._R8/(1024._R8*1024._R8)), ' MB '
 ALLOCATE( PMAT( NPROFS, MA, MA ) )
 PMAT = 0._R8
 ALLOCATE( RTMP( MA, MA ) )
 
 LOOPPROF : DO JPRF=1,NPROFS
    IF( NOBSPT( JPRF ) .EQ. 0 ) CYCLE LOOPPROF
    IF( NOBSPT( JPRF ) .EQ. 1 ) THEN
        PMAT(JPRF,1,1) = INS%ERR(JOSTART(JPRF))*INS%ERR(JOSTART(JPRF))
        CYCLE LOOPPROF
    ENDIF
    KO=0
    WRITE(1561+MYPROC,*) '///',JPRF,NOBSPT(JPRF),JOSTART(JPRF),JOEND(JPRF)
    LOOPINS2 : DO JO=JOSTART(JPRF),JOEND(JPRF)
       IF(INS%FLC(JO) .EQ. 0 ) CYCLE LOOPINS2
          KO=KO+1
          IF( KO .GT. MA ) CALL ABOR1('PREP_INSRINV: PMAT EXCEEDED 1')
          PMAT(JPRF,KO,KO) = INS%ERR(JO)*INS%ERR(JO)
          KO2=0
          LOOPINS3 : DO JO2=JOSTART(JPRF),JO-1

             IF(INS%FLC(JO2) .EQ. 0 ) CYCLE LOOPINS3
             KO2=KO2+1
             IF( KO2 .GT. MA ) CALL ABOR1('PREP_INSRINV: PMAT EXCEEDED 2')

             ! COMPUTE DISTANCE
             DIST=ABS(RODENS(JO)-RODENS(JO2))
             CR2 = RDSC(JPRF)

             ! COMPUTE CORRELATION
             IF(INS%PAR(JO).EQ.INS%PAR(JO2)) THEN
               COR = EXP(-0.5*(DIST*DIST)/(CR2*CR2))
             ELSE
               COR = 0._R8
             ENDIF

             ! COMPUTE COVARIANCE
             PMAT(JPRF,KO2,KO) = INS%ERR(JO)*INS%ERR(JO2)*COR
 
             ! TRUNCATE IF NEEDED
             IF(PMAT(JPRF,KO2,KO).LT.1.E-99_R8) PMAT(JPRF,KO2,KO) = 0._R8

             ! SYMMETRIZE
             PMAT(JPRF,KO,KO2) = PMAT(JPRF,KO2,KO)
          ENDDO LOOPINS3
    ENDDO LOOPINS2
 ENDDO LOOPPROF

 DEALLOCATE( RODENS)

! FOR THE TIME BEING, OMP PARALLELIZ. IS DISABLED
!...#ifdef SHARED_MEMORY
!...$OMP PARALLEL DEFAULT(SHARED), PRIVATE(JPRF,ZCHK,RTMP)
!...$OMP DO SCHEDULE(DYNAMIC,1)
!...#endif
 LOOPPROF2 : DO JPRF=1,NPROFS
    RTMP(:,:) = PMAT(JPRF,:,:)
    WRITE(1561+MYPROC,*) ' INVERSION: ',JPRF, NOBSPT(JPRF)
    IF ( NOBSPT(JPRF) .EQ. 1 ) THEN
       PMAT(JPRF,1,1) = 1._R8 / PMAT(JPRF,1,1)
    ELSE
       CALL MYSYMMINV(PMAT(JPRF,1:NOBSPT(JPRF),1:NOBSPT(JPRF)),NOBSPT(JPRF))
    ENDIF
    ZCHK = SUM ( MATMUL(RTMP(1:NOBSPT(JPRF),1:NOBSPT(JPRF)),PMAT(JPRF,1:NOBSPT(JPRF),1:NOBSPT(JPRF))) )
    WRITE(1561+MYPROC,*) ' RESULT   : ',JPRF, NOBSPT(JPRF),ZCHK, &
    & 100._R8*ABS(REAL(NOBSPT(JPRF),R8)-ZCHK)/REAL(NOBSPT(JPRF),R8)
    IF(LLDEB.AND.JPRF.EQ.1) THEN
       WRITE(CFMT,'(A,I3,A)') '(',NOBSPT(JPRF),'E15.6)'
       OPEN(938,FILE='dmat.dat')
       DO JO=1,NOBSPT(JPRF)
          WRITE(938,TRIM(CFMT)) RTMP(JO,1:NOBSPT(JPRF))
       ENDDO
       CLOSE(938)
       OPEN(939,FILE='imat.dat')
       DO JO=1,NOBSPT(JPRF)
          WRITE(939,TRIM(CFMT)) PMAT(JPRF,JO,1:NOBSPT(JPRF))
       ENDDO
       CLOSE(939)
       DO JO=JOSTART(JPRF),JOEND(JPRF)
          WRITE(939,'(I2)') INS%PAR(JO)
       ENDDO
    ENDIF
 ENDDO LOOPPROF2
!...#ifdef SHARED_MEMORY
!...$OMP END DO
!...$OMP END PARALLEL
!...#endif

 DEALLOCATE( RTMP )

 WRITE(IOUNLOG,*) ' NUMBER OF PROFS ',NPROFS
 WRITE(IOUNLOG,*) ' MIN/MAX NUMBER OF OBS PER PROF ', &
 & MINVAL( NOBSPT(1:NPROFS) ), MA

 CALL MYFRTPROF_WALL('PREP_INSRINV: PREPARING FOR INS R INVERSION',1)

END SUBROUTINE PREP_INSRINV

SUBROUTINE MYSYMMINV(A,N)

! INVERSE OF SYMMETRIC MATRIX VIA CHOLESKY DECOMP.
! USE LAPACK STANDARD ROUTINES
! A.S. - 2016

  IMPLICIT NONE
  INTEGER(I4), INTENT(IN) :: N
  REAL(R8), INTENT(INOUT) :: A(N,N)
  REAL(R8) :: A2(N,N)
  INTEGER(I4) :: IERR,I,J
  IF ( N .LE. 0 ) THEN
     CALL ABOR1('SYMMINV: 0-LENGTH ARRAY')
  ENDIF
  IF ( N .EQ. 1 ) THEN
     A(1,1) = 1._R8 / A(1,1)
     RETURN
  ENDIF
  A2=A
  IF (R8 == KIND(1._8)) THEN
     CALL DPOTRF('L',N,A,N,IERR)
     IF (IERR .NE. 0) THEN
        WRITE(IOUNLOG,*) ' DPOTRF RETURNED ', IERR
        WRITE(IOUNLOG,*) ' FAILED INVERSION WITH CHOLESKY DECOMP., DIM = ',N
        WRITE(IOUNLOG,*) ' USING EIGEN-DECOMP. PSEUDO-INVERSE, DIM = ',N
        A=A2
        CALL MYPINV(N,A)
        RETURN
     END IF
     CALL DPOTRI('L',N,A,N,IERR)
     IF (IERR .NE. 0) THEN
        WRITE(IOUNERR,*) ' DPOTRI RETURNED ', IERR
        WRITE(IOUNLOG,*) ' FAILED INVERSION, DIM = ',N
        DO J=1,N
         DO I=1,N
          WRITE(IOUNLOG,*) I, J, A(I,J)
         ENDDO
        ENDDO
        CALL ABOR1('SYMMINV: ERROR IN DPOTRI')
     END IF
  ELSE IF (R8 == KIND(1._4)) THEN
     CALL SPOTRF('L',N,A,N,IERR)
     IF (IERR .NE. 0) THEN
        WRITE(IOUNERR,*) ' SPOTRF RETURNED ', IERR
        WRITE(IOUNLOG,*) ' FAILED INVERSION, DIM = ',N
        DO J=1,N
         DO I=1,N
          WRITE(IOUNLOG,*) I, J, A(I,J)
         ENDDO
        ENDDO
        CALL ABOR1('SYMMINV: ERROR IN SPOTRF')
     END IF
     CALL SPOTRI('L',N,A,N,IERR)
     IF (IERR .NE. 0) THEN
        WRITE(IOUNERR,*) ' SPOTRI RETURNED ', IERR
        WRITE(IOUNLOG,*) ' FAILED INVERSION, DIM = ',N
        DO J=1,N
         DO I=1,N
          WRITE(IOUNLOG,*) I, J, A(I,J)
         ENDDO
        ENDDO
        CALL ABOR1('SYMMINV: ERROR IN SPOTRI')
     END IF
  ELSE
     CALL ABOR1('SYMMINV: UNRECOGNIZED KIND')
  ENDIF
  DO J=2,N
   DO I=1,J-1
     A(I,J) = A(J,I)
   ENDDO
  ENDDO
END SUBROUTINE MYSYMMINV

subroutine pinv2 (N, A)

      integer(i4), intent(in) :: N
      real(R8), intent (inout) :: A(N,N)
      real(R8), parameter      :: tolerance = 1.E-6_R8
      real(R8), allocatable    :: work(:)
      real(R8) :: B(N,N)
      integer(i4) :: IPIV(N),j,INFO,LW,I

      B = 0._R8
      DO J =1,N
         B(J,J) = 1._R8
      ENDDO

      LW = N*N
      ALLOCATE(WORK(LW))
      CALL DGESV(N,N,A,N,IPIV,B,N,INFO)
      DEALLOCATE(WORK)

      WRITE(1561+MYPROC,*) 'DGESV RETURNED ',INFO

      A = B
      !DO J=2,N
      ! DO I=1,J-1
      !   A(I,J) = A(J,I)
      ! ENDDO
      !ENDDO

END subroutine pinv2

function pinv_workspace(N) result(lwork)
!     Determine the workspace needed for dsyev.
      implicit none
      integer(i4),intent(in) :: N
      integer(i4) :: lwork
      integer(i4) :: info
      real(R8) :: dummy,rwork
      call dsyev('V','U', N, dummy,N, dummy, rwork, -1, info )
      lwork = ceiling(rwork)
end function pinv_workspace


subroutine mypinv (N, A)

!     Compute pseudo-inverse A+ of symmetric A in factored form U D+ U', where
!     U overwrites A and D is diagonal matrix D+.

!     Saunders notes: Working with the factors of the pseudo-inverse is
!                     preferable to multiplying them together (more stable,
!                     less arithmetic).  Also, the choice of tolerance is
!                     not straightforward.  If A is noisy, try 1.e-2 * || A ||,
!                     else machine-eps. * || A ||  (1 norms).
!     Arguments:

      integer(i4), intent(in) :: N
      real(R8), intent (inout) :: A(N,N)
      real(R8), parameter      :: tolerance = 1.490116e-11_R8
      real(R8), allocatable    :: work(:), D(:)
      real(R8)  :: C(N,N)

!     Local variables:

      integer :: i, info, NWS, J, K

!     Execution:

      NWS = PINV_WORKSPACE(N) 

      ALLOCATE( work(NWS) )
      ALLOCATE( D(N) )

!     Eigendecomposition/SVD of symmetric A:

      call dsyev ('V', 'U', N, A, N, D, work, NWS, info)

!     Diagonal factor D+ of pseudo-inverse:

      do i = 1, N
         if (D(i) > tolerance) then
            D(i) = 1._R8 / D(i)
         else
            D(i) = 0._R8
         end if
      end do

      DO K=1,N
        DO J=1,N
          C(J,K) = 0._R8
          DO I=1,N
            C(J,K) = C(J,K) + A(J,I) * D(I) * A(K,I)
          END DO
        END DO
      END DO

      DEALLOCATE( work, D)

      A = C

      end subroutine mypinv

SUBROUTINE READ_DENS(I,J,K,TMP)
  IMPLICIT NONE
  INTEGER(I4), INTENT(IN) :: I,J,K
  REAL(R8), INTENT(OUT) :: TMP(I,J,K)
  REAL(R8) :: TMPT(I,J,K), TMPS(I,J,K)
  INTEGER(I4) :: JI,JJ,JK
  CALL GETNCVAR('DENSITY.nc','votemper',GRD%IM,GRD%JM,GRD%KM,TMPT)
  CALL GETNCVAR('DENSITY.nc','vosaline',GRD%IM,GRD%JM,GRD%KM,TMPS)
  DO JK=1,K
   DO JJ=1,J
    DO JI=1,I
      TMP(JI,JJ,JK) = RHO_UNESCO(MAX(TMPS(JI,JJ,JK),0._R8),TMPT(JI,JJ,JK),0._R8,.FALSE.)
    ENDDO
   ENDDO
  ENDDO
  WHERE( .NOT. ABS(TMP) .LT. 1.E+9 ) TMP = 0._R8
  CALL WRITE_VARNCDF('DENSITY_'//TRIM(CMPIDOM)//'.nc',I,J,K,TMP,&
  & 'vodensit','kg/m3','Seawater density')
END SUBROUTINE READ_DENS

END MODULE PROF_CORREL
