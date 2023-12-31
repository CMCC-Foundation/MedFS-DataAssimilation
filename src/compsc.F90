SUBROUTINE COMPSC
USE RECFILTER

!-- COMPUTE SCALING FACTOR FOR RF

USE GRD_STR

IMPLICIT  NONE

INTEGER(KIND=I4) :: NX,NY,NZ
INTEGER(KIND=I4) :: NSPL,NTB
REAL(KIND=R8) :: DSMN,DSMX,DSL,DST,E
INTEGER(KIND=I4) :: JLEV,K,J,I
REAL(KIND=R8),ALLOCATABLE :: SFCT(:),AL(:),BT(:)
INTEGER(KIND=I4),ALLOCATABLE :: JNXX(:),IAUX(:)
REAL(KIND=R8),ALLOCATABLE :: ALNS(:,:),SCNS(:,:)

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) ' /// COMPUTE RF SCALING FACTORS'
WRITE(IOUNLOG,*)

#ifndef USE_POINTERS
IF(.NOT.ALLOCATED(GRD%DX).OR..NOT.ALLOCATED(GRD%DY)) THEN
   WRITE(IOUNERR,*) 'COMPSC CALLED BEFORE GRID INITIALIZATION'
   CALL ABOR1('ERROR IN COMPSC, CANNOT PROCEED')
ENDIF
#endif

IF(.NOT.LLRFINIT) &
& CALL ABOR1('COMPSC CALLED BEFORE RF INITIALIZATION')

NX=GRD%IM
NY=GRD%JM
NZ=GRD%KM

ALLOCATE(SCXNS(NX,NY,NTOTV), SCYNS(NX,NY,NTOTV))
ALLOCATE(IAUX(1))

NSPL=MAX(NX,NY)
ALLOCATE (SFCT(NSPL),JNXX(NSPL),AL(NSPL),BT(NSPL) )
NTB = MIN(20,MIN(NX,NY))

ALLOCATE(ALNS(NTB,NTOTV))
ALLOCATE(SCNS(NTB,NTOTV))

DSMN=MIN(MINVAL(GRD%DX),MINVAL(GRD%DY))
DSMX=MAX(MAXVAL(GRD%DX),MAXVAL(GRD%DY))

WRITE(IOUNLOG,*) ' MAX,MIN DIM.:',NSPL,NTB
WRITE(IOUNLOG,*) ' MAX,MIN RES.:',DSMX,DSMN

DSMX=DSMX+MAX(1._R8,(DSMX-DSMN)/(REAL(NTB,KIND=R8)-2._R8))
DSL=(DSMX-DSMN)/(REAL(NTB,KIND=R8)-1._R8)

! (NTOTV EVENTUALLY ACCOUNTS FOR BOTH T,S)
DO JLEV=1,NTOTV
  DO K=1,NTB
    DST=DSMN+(K-1._R8)*DSL
    E=(2._R8*RF_NTR) * DST**2 / (4._R8 * ZGA_CRAD(JLEV)**2)
    ALNS(K,JLEV) = 1._R8 + E - SQRT(E*(E+2._R8))
    SFCT(:)=0._R8
    AL(:)=ALNS(K,JLEV)
    BT(:)=ALNS(K,JLEV)
    DO J=1,NSPL
      JNXX(J)=J
    ENDDO
    SFCT(NSPL/2 + 1)=1._R8
!... NOTE THAT USE OF IAUX AVOIDS XLF ERRORS AS IT SHOULD
!    BE A VECTOR AND NOT A SCALAR
    IAUX(1) = NSPL
    CALL RCFL_Y_OLD(1,NSPL,1,NSPL,AL,BT,SFCT,JNXX,IAUX)
    CALL RCFL_Y_AD_OLD(1,NSPL,1,NSPL,AL,BT,SFCT,JNXX,IAUX)
    SCNS(K,JLEV) = SFCT(NSPL/2+1)
  ENDDO
ENDDO
DEALLOCATE (SFCT,JNXX,AL,BT)

DO JLEV=1,NTOTV
   DO J=1,NY
     DO I=1,NX
         DST = ( GRD%DX(I,J) - DSMN )/DSL
         K = INT(DST) + 1
         DST = DST - REAL(K-1)
         SCXNS(I,J,JLEV) = &
         & SQRT( 1._R8/ (SCNS(K,JLEV)*(1._R8-DST) + SCNS(K+1,JLEV)*DST) )
         DST = ( GRD%DY(I,J) - DSMN )/DSL
         K = INT(DST) + 1
         DST = DST - REAL(K-1)
         SCYNS(I,J,JLEV) = &
         & SQRT( 1._R8/ (SCNS(K,JLEV)*(1._R8-DST) + SCNS(K+1,JLEV)*DST) )
     ENDDO
   ENDDO
   IF (LL_PRINT_RECFILT) THEN
      WRITE(IOUNLOG,*) 'SCXNS ', JLEV, ': ',SUM(SCXNS(:,:,JLEV))/(NY*NX),&
      & MINVAL(SCXNS(:,:,JLEV)),MAXVAL(SCXNS(:,:,JLEV))
      WRITE(IOUNLOG,*) 'SCYNS ', JLEV, ': ',SUM(SCYNS(:,:,JLEV))/(NY*NX),&
      & MINVAL(SCYNS(:,:,JLEV)),MAXVAL(SCYNS(:,:,JLEV))
   ENDIF
ENDDO

DEALLOCATE(ALNS,SCNS)

RETURN
END SUBROUTINE COMPSC
