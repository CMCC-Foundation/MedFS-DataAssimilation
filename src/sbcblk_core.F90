MODULE BLKOCE_CORE

 USE TLAD_VARS
 USE MYFRTPROF

 IMPLICIT NONE
 REAL(WP), DIMENSION(:,:), ALLOCATABLE :: TEM0, WND, SAL0, CE, CH, Q10, EMP

CONTAINS

   SUBROUTINE BLK_OCE_CORE(JPI,JPJ,TEM0_TL,SAL0_TL,TEM0A,SAL0A)
      INTEGER, INTENT(IN) :: JPI, JPJ
      REAL(WP) , INTENT(IN), DIMENSION(JPI,JPJ) :: TEM0_TL,SAL0_TL
      REAL(WP) , INTENT(INOUT), DIMENSION(JPI,JPJ) :: TEM0A, SAL0A
      !
      INTEGER  ::   JI, JJ               ! DUMMY LOOP INDICES
      REAL(WP), DIMENSION(JPI,JPJ) :: ZQLW_TL,ZQSATW_TL,ZAUX,ZQSB_TL,ZQLA_TL,ZEVAP_TL
      REAL(WP), DIMENSION(JPI,JPJ) :: QNS_TL,EMP_TL,ZQSATW
      CALL MYFRTPROF_WALL('BLK_OCE_CORE: TL BULK FLX',0)
      ! ----------------------------------------------------------------------------- !
      !      I   RADIATIVE FLUXES                                                     !
      ! ----------------------------------------------------------------------------- !
      ZQLW_TL(:,:) = - (4._WP*STEF*TEM0(:,:)*TEM0(:,:)*TEM0(:,:)*TEM0_TL(:,:))*TMASK(:,:,1)
      ! ----------------------------------------------------------------------------- !
      !     II    TURBULENT FLUXES                                                    !
      ! ----------------------------------------------------------------------------- !
      ZQSATW(:,:) = ZCOEF_QSATW*EXP(-5107.4_R8/TEM0(:,:))
      ZQSATW_TL(:,:) = ZCOEF_QSATW*5107.4_R8*TEM0_TL(:,:)*EXP(-5107.4_R8/TEM0(:,:))/TEM0(:,:)**2
      !
      ZEVAP_TL = 0._R8
      ZAUX   = RHOA  *CE(:,:)* ( ZQSATW(:,:) - Q10(:,:) ) * WND(:,:) 
      WHERE ( ZAUX .GT. 0._R8 ) ZEVAP_TL = RHOA *CE(:,:) *  ZQSATW_TL(:,:) * WND(:,:) 
      ZQSB_TL (:,:) =  RHOA*CPA*CH(:,:)*TEM0_TL(:,:) * WND(:,:)     ! SENSIBLE HEAT
      ZQLA_TL (:,:) = LV * ZEVAP_TL(:,:)         ! LATENT HEAT
      ! ----------------------------------------------------------------------------- !
      !     III    TOTAL FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
      QNS_TL(:,:) = ZQLW_TL(:,:) - ZQSB_TL(:,:) - ZQLA_TL(:,:)      ! DOWNWARD NON SOLAR FLUX
      EMP_TL(:,:) = ZEVAP_TL(:,:) 
      ! ----------------------------------------------------------------------------- !
      !      IV    TENDENCIES                                                         !
      ! ----------------------------------------------------------------------------- !
      TEM0A = TEM0A + TMASK(:,:,1)* RO0CPR * QNS_TL / E3T(1)
      SAL0A = SAL0A + TMASK(:,:,1)* ZSRAU * (EMP_TL*SAL0 + EMP*SAL0_TL) / E3T(1)
      CALL MYFRTPROF_WALL('BLK_OCE_CORE: TL BULK FLX',1)
      !
   END SUBROUTINE BLK_OCE_CORE

   SUBROUTINE BLK_OCE_CORE_AD(JPI,JPJ,TEM0_AD,SAL0_AD,TEM0A,SAL0A)
      INTEGER, INTENT(IN) :: JPI, JPJ
      REAL(WP) , INTENT(INOUT), DIMENSION(JPI,JPJ) :: TEM0_AD,SAL0_AD
      REAL(WP) , INTENT(IN), DIMENSION(JPI,JPJ) :: TEM0A, SAL0A
      !
      INTEGER  ::   JI, JJ               ! DUMMY LOOP INDICES
      REAL(WP), DIMENSION(JPI,JPJ) :: EMP_AD, QNS_AD, EVAP_AD, QLW_AD, QSB_AD, QLA_AD
      REAL(WP), DIMENSION(JPI,JPJ) :: ZQSATW, ZAUX, ZQSATW_AD

      CALL MYFRTPROF_WALL('BLK_OCE_CORE_AD: AD BULK FLX',0)

      ! ----------------------------------------------------------------------------- !
      !      IV    TENDENCIES                                                         !
      ! ----------------------------------------------------------------------------- !
      EMP_AD =  TMASK(:,:,1)* ZSRAU * SAL0 * SAL0A / E3T(1)
      SAL0_AD = SAL0_AD + TMASK(:,:,1)* ZSRAU * EMP * SAL0A / E3T(1)
      QNS_AD =  TMASK(:,:,1)* RO0CPR * TEM0A / E3T(1)
      ! ----------------------------------------------------------------------------- !
      !     III    TOTAL FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
      EVAP_AD= EMP_AD
      QLW_AD = QNS_AD
      QSB_AD = -QNS_AD
      QLA_AD = -QNS_AD
      ! ----------------------------------------------------------------------------- !
      !     II    TURBULENT FLUXES                                                    !
      ! ----------------------------------------------------------------------------- !
      TEM0_AD = TEM0_AD + QSB_AD * RHOA*CPA*CH(:,:)*WND(:,:)
      EVAP_AD = EVAP_AD + LV*QLA_AD

      ZQSATW(:,:) = ZCOEF_QSATW*EXP(-5107.4_R8/TEM0(:,:))
      ZAUX     = RHOA    *CE(:,:)* ( ZQSATW(:,:) - Q10(:,:) ) * WND(:,:) 
      ZQSATW_AD = 0._R8
      WHERE ( ZAUX .GT. 0._R8 ) ZQSATW_AD = RHOA*CE(:,:) *EVAP_AD*WND(:,:) 
      TEM0_AD = TEM0_AD + EXP(-(5107.4_R8/TEM0))*ZCOEF_QSATW*5107.4_R8*ZQSATW_AD/TEM0**2
      ! ----------------------------------------------------------------------------- !
      !      I   RADIATIVE FLUXES                                                     !
      ! ----------------------------------------------------------------------------- !
      TEM0_AD = TEM0_AD - &
      & (4._WP*STEF*TEM0(:,:)*TEM0(:,:)*TEM0(:,:)*QLW_AD(:,:))*TMASK(:,:,1)

      CALL MYFRTPROF_WALL('BLK_OCE_CORE_AD: AD BULK FLX',1)

   END SUBROUTINE BLK_OCE_CORE_AD

END MODULE BLKOCE_CORE
