SUBROUTINE SRFBAL_AD(RDT,TEMB,TEM_AD,CE,CH,BLH,Q10,WND,T_AML_AD,Q_AML_AD)

USE SET_KND
USE GRD_STR

IMPLICIT NONE

REAL(R8), INTENT(IN)  :: RDT
REAL(R8), DIMENSION(GRD%IM,GRD%JM), INTENT(IN)  :: TEMB
REAL(R8), DIMENSION(GRD%IM,GRD%JM), INTENT(INOUT)  :: TEM_AD
REAL(R8), DIMENSION(GRD%IM,GRD%JM), INTENT(IN)  :: CE, CH, BLH, Q10, WND
REAL(R8), DIMENSION(GRD%IM,GRD%JM), INTENT(IN) :: T_AML_AD, Q_AML_AD

REAL(R8), ALLOCATABLE, DIMENSION(:,:) :: ZQLW_AD, ZQSATW, ZQSATW_AD
REAL(R8), ALLOCATABLE, DIMENSION(:,:) :: ZEVA_AD, ZAUX
REAL(R8), ALLOCATABLE, DIMENSION(:,:) :: ZQSB_AD, ZQLA_AD, TEM0

REAL(R8), PARAMETER :: STEF =5.67E-8_R8
REAL(R8), PARAMETER :: RHOA =    1.22_R8
REAL(R8), PARAMETER :: CPA  = 1000.5_R8
REAL(R8), PARAMETER :: LV   =    2.5E6_R8
REAL(R8), PARAMETER :: RT0  =  273.15_R8
REAL(R8), PARAMETER :: ZCOEF_QSATW = 0.98_R8 * 640380._R8 / RHOA

ALLOCATE( ZQLW_AD(GRD%IM,GRD%JM) )
ALLOCATE( ZQSATW (GRD%IM,GRD%JM) )
ALLOCATE( ZQSATW_AD (GRD%IM,GRD%JM) )
ALLOCATE( ZEVA_AD (GRD%IM,GRD%JM) )
ALLOCATE( ZAUX   (GRD%IM,GRD%JM) )
ALLOCATE( ZQSB_AD(GRD%IM,GRD%JM) )
ALLOCATE( ZQLA_AD(GRD%IM,GRD%JM) )
ALLOCATE( TEM0   (GRD%IM,GRD%JM) )

TEM0 = TEMB + RT0

! 2M PARAMETERS
ZQLW_AD    = ( T_AML_AD * RDT / (RHOA*CPA*BLH) ) * GRD%MSK(:,:,1)
ZQSB_AD    = ( T_AML_AD * RDT / (RHOA*CPA*BLH) ) * GRD%MSK(:,:,1)
ZEVA_AD    = ( Q_AML_AD * RDT / (RHOA*BLH)     ) * GRD%MSK(:,:,1)

! SENSIBLE HEAT
TEM_AD = TEM_AD + ZQSB_AD*RHOA*CPA*CH*WND

!EVAPORATION
ZQSATW_AD = 0._R8
ZAUX   = RHOA  *CE* ( ZQSATW - Q10 ) * WND
WHERE ( ZAUX .GT. 0._R8 ) ZQSATW_AD = ZEVA_AD * RHOA *CE * WND
TEM_AD = TEM_AD + ZQSATW_AD*ZCOEF_QSATW*5107.4_R8*EXP(-5107.4_R8/TEM0)/TEM0**2

! UPWELLING LONGWAVE
TEM_AD = TEM_AD + (4._R8*STEF*TEM0*TEM0*TEM0*ZQLW_AD)*GRD%MSK(:,:,1)

DEALLOCATE( ZQLW_AD, ZQSATW, ZQSATW_AD, ZEVA_AD, ZAUX, ZQSB_AD, ZQLA_AD, TEM0)

END SUBROUTINE SRFBAL_AD
