SUBROUTINE COMPALPHA_LV

USE SET_KND
USE RECFILTER
USE GRD_STR
USE LBCLNK
USE MYFRTPROF
USE MYNETCDF

!-- COMPUTE ALPHA COEFFICIENTS FOR RF

IMPLICIT  NONE

INTEGER(KIND=I4) :: NX,NY,NZ
INTEGER(KIND=I4) :: I,J,JLEV,IERR
REAL(KIND=R8) :: E,DST, TCR, TCRM, SIGMA

!INTEGER(KIND=I4) :: ITMP(GRD%IM,GRD%JM)

CALL MYFRTPROF_WALL('COMPALPHA_LV: COMPUTE ALPHA COEFFS',0)

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) ' /// COMPUTE ALPHA COEFFICIENTS'
WRITE(IOUNLOG,*)

CALL FLUSH(IOUNLOG)

#ifndef USE_POINTERS
IF(.NOT.ALLOCATED(NEWDX).OR..NOT.ALLOCATED(NEWDY)) THEN
   WRITE(IOUNERR,*) 'COMPALPHA CALLED BEFORE GRID INITIALIZATION'
   CALL ABOR1('ERROR IN COMPALPHA, CANNOT PROCEED')
ENDIF
#endif

IF(.NOT.LLRFINIT) &
& CALL ABOR1('COMPALPHA CALLED BEFORE RF INITIALIZATION')

NX=GRD%IM
NY=GRD%JM
NZ=GRD%KM

WRITE(IOUNLOG,*) 'DIMENSIONS ARE ', NX, NY, NTOTV
CALL FLUSH(IOUNLOG)

ALLOCATE(ALXNS(NX,NY,NTOTV),&
       & BTXNS(NX,NY,NTOTV),&
       & ALYNS(NX,NY,NTOTV),&
       & BTYNS(NX,NY,NTOTV),STAT=IERR)

ALXNS = 0._R8
BTXNS = 0._R8
ALYNS = 0._R8
BTYNS = 0._R8

IF( RCF_TYPE .EQ. 2 ) THEN
   ALLOCATE(GMXNS(NX,NY,NTOTV),&
          & DLXNS(NX,NY,NTOTV),&
          & GMYNS(NX,NY,NTOTV),&
          & DLYNS(NX,NY,NTOTV),STAT=IERR)
   GMXNS = 0._R8
   DLXNS = 0._R8
   GMYNS = 0._R8
   DLYNS = 0._R8
ENDIF

IF(IERR .NE. 0 ) THEN
   CALL ABOR1('ALLOCATION PROBLEMS IN COMPALPHA_LV')
ENDIF

WRITE(IOUNLOG,*) 'ALLOCATION DONE'
CALL FLUSH(IOUNLOG)

DO JLEV=1,NTOTV

 IF (LL_PRINT_RECFILT) THEN
  WRITE(IOUNLOG,*) ' DOING LEVEL ', JLEV, ' OVER ', NTOTV
  CALL FLUSH(IOUNLOG)
 ENDIF

 IF( RCF_TYPE .EQ. 1 ) THEN

   DO J=1,NY
      DO I=2,NX

         TCR = CRX(I,J,JLEV)
         TCRM= CRX(I-1,J,JLEV)

         DST = (NEWDX(I-1,J) + NEWDX(I,J)) * 0.5_R8
         TCR = (TCR + TCRM) * 0.5_R8
         E   = (2._R8 * RF_NTR) * DST**2 / (4._R8 *TCR**2)
         ALXNS(I,J,JLEV) = 1._R8 + E - SQRT(E*(E+2._R8))

      ENDDO

      DO I=1,NX-1

         TCR = CRX(I,J,JLEV)
         TCRM= CRX(I+1,J,JLEV)

         TCR = (TCR + TCRM) * 0.5_R8
         DST = (NEWDX(I,J) + NEWDX(I+1,J)) * 0.5_R8
         E   = (2._R8 * RF_NTR) * DST**2 / (4._R8 *TCR**2)
         BTXNS(I,J,JLEV) = 1._R8 + E - SQRT(E*(E+2._R8))

      ENDDO

   ENDDO

   DO J=2,NY
      DO I=1,NX

         TCR = CRY(I,J,JLEV)
         TCRM= CRY(I,J-1,JLEV)

         TCR = (TCR + TCRM) * 0.5_R8
         DST = (NEWDY(I,J-1) + NEWDY(I,J)) * 0.5_R8
         E   = (2._R8 * RF_NTR) * DST**2 / (4._R8 *TCR**2)
         ALYNS(I,J,JLEV) = 1._R8 + E - SQRT(E*(E+2._R8))

      ENDDO
   ENDDO

   DO J=1,NY-1
      DO I=1,NX

         TCR = CRY(I,J,JLEV)
         TCRM= CRY(I,J+1,JLEV)

         TCR = (TCR + TCRM) * 0.5_R8
         DST = (NEWDY(I,J) + NEWDY(I,J+1)) * 0.5_R8
         E   = (2._R8 * RF_NTR) * DST**2 / (4._R8 *TCR**2)
         BTYNS(I,J,JLEV) = 1._R8 + E - SQRT(E*(E+2._R8))

      ENDDO
   ENDDO

 ELSEIF( RCF_TYPE .EQ. 2 ) THEN

   DO J=1,NY
      DO I=1,NX
         TCR = CRX(I,J,JLEV)
         SIGMA = MAX(0.5_R8,TCR/GRD%DX(I,J) )
         CALL DEF_COEF(SIGMA,ALXNS(I,J,JLEV),BTXNS(I,J,JLEV),&
         & GMXNS(I,J,JLEV),DLXNS(I,J,JLEV))

         TCR = CRY(I,J,JLEV)
         SIGMA = MAX(0.5_R8,TCR/GRD%DY(I,J) )
         CALL DEF_COEF(SIGMA,ALYNS(I,J,JLEV),BTYNS(I,J,JLEV),&
         & GMYNS(I,J,JLEV),DLYNS(I,J,JLEV))
      ENDDO
   ENDDO

 ENDIF

ENDDO

!CALL WRITE_VARNCDF('alx.nc',GRD%IM,GRD%JM,ALXNS(:,:,1),'alxns','','')
!CALL WRITE_VARNCDF('aly.nc',GRD%IM,GRD%JM,ALYNS(:,:,1),'alyns','','')
!CALL WRITE_VARNCDF('btx.nc',GRD%IM,GRD%JM,BTXNS(:,:,1),'btxns','','')
!CALL WRITE_VARNCDF('bty.nc',GRD%IM,GRD%JM,BTYNS(:,:,1),'btyns','','')
!ITMP = LLMSR(:,:,1)
!CALL WRITE_VARNCDF('llmsrs.nc',GRD%IM,GRD%JM,ITMP,'lsm','','')
!IF( RCF_TYPE .EQ. 2 ) THEN
! CALL WRITE_VARNCDF('gmx.nc',GRD%IM,GRD%JM,GMXNS(:,:,1),'gmxns','','')
! CALL WRITE_VARNCDF('gmy.nc',GRD%IM,GRD%JM,GMYNS(:,:,1),'gmyns','','')
! CALL WRITE_VARNCDF('dlx.nc',GRD%IM,GRD%JM,DLXNS(:,:,1),'dlxns','','')
! CALL WRITE_VARNCDF('dly.nc',GRD%IM,GRD%JM,DLYNS(:,:,1),'dlyns','','')
!ENDIF

CALL MYFRTPROF_WALL('COMPALPHA_LV: COMPUTE ALPHA COEFFS',1)

END SUBROUTINE COMPALPHA_LV
