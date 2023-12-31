MODULE TLAD_TEST

USE SET_KND
USE GRD_STR
USE TLAD_VARS
USE TRAADV_CEN2H
USE TRAADV_CEN2V
USE TRALDF_LAP
USE TRAZDF_EXP
USE TRAZDF_IMP
USE IOUNITS
USE MYFRTPROF
USE TLAD
USE BLKOCE_CORE
USE ADTEST
USE WEAKLY , ONLY : LL_WEAKLY

IMPLICIT NONE

CONTAINS

SUBROUTINE TLADMOD_TEST

IMPLICIT NONE

REAL(R8), DIMENSION(JPI,JPJ,JPK) :: GAT, GAS, GTEM, GTEM_AD, GAT0,GAS0
REAL(R8), DIMENSION(JPI,JPJ) :: GTEM2_AD, GSAL2_AD, GTEM2, GSAL2, &
& GAT2, GAS2, GAT2B, GAS2B, GAT2B0, GAS2B0
REAL(R8), DIMENSION(JPI,JPJ,JPK) :: TBK, SBK
REAL(R8), ALLOCATABLE,DIMENSION(:) :: ZV1,ZV2,ZV3,ZV4
REAL(R8) :: ZT1, ZT2
INTEGER(I4) :: JT, JO, JK, JI, JJ, KK
LOGICAL :: LLSAVE

 CALL MYFRTPROF_WALL('TLMOD_TEST: TEST ADJOINT',0)

GAT = 0._R8
GAS = 0._R8

   JT=1

   TLAD_TS = JT

   WRITE(IOUNLOG,*) ' TL TEST TIMESTEP : ',JT

   CALL SET_BGF1(JT)
   CALL SET_BGF2(JT)
   CALL SET_BGF3(JT)
   WRITE(IOUNLOG,*) ' MIN/MAX UBT ',MINVAL(UBT),MAXVAL(UBT),MINVAL(WBT),MAXVAL(WBT)

   !... TRACER ADVECTION
   CALL RANDOM_FLD(GTEM)
   CALL RANDOM_FLD(GAS)
   GAT = 0._R8
   GTEM_AD = 0._R8
   CALL TRA_ADV_CEN2H(JPI,JPJ,JPK,UBT,VBT,GTEM,GAT)
   CALL TRA_ADV_CEN2H_AD(JPI,JPJ,JPK,UBT,VBT,GTEM_AD,GAS)
   ZT1 = MDOT_PRODUCT(GTEM,GTEM_AD)
   ZT2 = MDOT_PRODUCT(GAT,GAS)
   WRITE(IOUNLOG,*) ' ADJOINT TEST FOR TRA_ADV_CEN2H'
   WRITE(IOUNLOG,*) ' MIN/MAX UBT ',MINVAL(UBT),MAXVAL(UBT),MINVAL(WBT),MAXVAL(WBT)
   WRITE(IOUNLOG,*) ' GTEM EXTRMS = ',MINVAL(GTEM),MAXVAL(GTEM),&
   &SUM(GTEM)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAS EXTRMS = ',MINVAL(GAS),MAXVAL(GAS),&
   &SUM(GAS)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GTEM_AD EXTRMS = ',MINVAL(GTEM_AD),MAXVAL(GTEM_AD),&
   &SUM(GTEM_AD)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAT EXTRMS = ',MINVAL(GAT),MAXVAL(GAT),&
   &SUM(GAT)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
   WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
   WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
   WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
   WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
   WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
   WRITE(IOUNLOG,*) 
   WRITE(IOUNLOG,*) 

   CALL RANDOM_FLD(GTEM)
   CALL RANDOM_FLD(GAS)
   GAT = 0._R8
   GTEM_AD = 0._R8
   CALL TRA_ADV_CEN2V(JPI,JPJ,JPK,WBT,GTEM,GAT)
   CALL TRA_ADV_CEN2V_AD(JPI,JPJ,JPK,WBT,GTEM_AD,GAS)
   ZT1 = MDOT_PRODUCT(GTEM,GTEM_AD)
   ZT2 = MDOT_PRODUCT(GAT,GAS)
   WRITE(IOUNLOG,*) ' ADJOINT TEST FOR TRA_ADV_CEN2V'
   WRITE(IOUNLOG,*) ' MIN/MAX UBT ',MINVAL(UBT),MAXVAL(UBT),MINVAL(WBT),MAXVAL(WBT)
   WRITE(IOUNLOG,*) ' GTEM EXTRMS = ',MINVAL(GTEM),MAXVAL(GTEM),&
   &SUM(GTEM)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAS EXTRMS = ',MINVAL(GAS),MAXVAL(GAS),&
   &SUM(GAS)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GTEM_AD EXTRMS = ',MINVAL(GTEM_AD),MAXVAL(GTEM_AD),&
   &SUM(GTEM_AD)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAT EXTRMS = ',MINVAL(GAT),MAXVAL(GAT),&
   &SUM(GAT)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
   WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
   WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
   WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
   WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
   WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
   WRITE(IOUNLOG,*)
   WRITE(IOUNLOG,*)

   !... TRACER DIFFUSION
   CALL RANDOM_FLD(GTEM)
   CALL RANDOM_FLD(GAS)
   GAT = 0._R8
   GTEM_AD = 0._R8
   CALL TRA_LDF_LAP(JPI,JPJ,JPK,GTEM,GAT)
   CALL TRA_LDF_LAP_AD(JPI,JPJ,JPK,GTEM_AD,GAS)
   ZT1 = MDOT_PRODUCT(GTEM,GTEM_AD)
   ZT2 = MDOT_PRODUCT(GAT,GAS)
   WRITE(IOUNLOG,*) ' ADJOINT TEST FOR TRA_LDF_LAP'
   WRITE(IOUNLOG,*) ' GTEM EXTRMS = ',MINVAL(GTEM),MAXVAL(GTEM),&
   &SUM(GTEM)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAS EXTRMS = ',MINVAL(GAS),MAXVAL(GAS),&
   &SUM(GAS)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GTEM_AD EXTRMS = ',MINVAL(GTEM_AD),MAXVAL(GTEM_AD),&
   &SUM(GTEM_AD)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAT EXTRMS = ',MINVAL(GAT),MAXVAL(GAT),&
   &SUM(GAT)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
   WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
   WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
   WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
   WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
   WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
   WRITE(IOUNLOG,*) 
   WRITE(IOUNLOG,*) 

   !... IMPLICIT VERTICAL DIFFUSION
   CALL RANDOM_FLD(GTEM)
!   CALL RANDOM_FLD(GAS,.TRUE.)
   CALL RANDOM_FLD(GAS)
   GAT = 0._R8
   GTEM_AD = 0._R8
   GAS0=GAS
   CALL TRA_ZDF_IMP(2*RDT,GTEM,GAT,AVBT)
   CALL TRA_ZDF_IMP_AD(2*RDT,GTEM_AD,GAS,AVBT)
   ZT1 = MDOT_PRODUCT(GTEM,GTEM_AD)
!   DO JK=1,JPK
!     GAS(:,:,JK)=GAS(:,:,JK)/(E1T*E2T*E3T(JK))
!   ENDDO
   ZT2 = MDOT_PRODUCT(GAT,GAS0)
   WRITE(IOUNLOG,*) ' ADJOINT TEST FOR TRA_ZDF_IMP'
   WRITE(IOUNLOG,*) ' GTEM EXTRMS = ',MINVAL(GTEM),MAXVAL(GTEM),&
   &SUM(GTEM)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAS EXTRMS = ',MINVAL(GAS),MAXVAL(GAS),&
   &SUM(GAS)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GTEM_AD EXTRMS = ',MINVAL(GTEM_AD),MAXVAL(GTEM_AD),&
   &SUM(GTEM_AD)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAT EXTRMS = ',MINVAL(GAT),MAXVAL(GAT),&
   &SUM(GAT)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
   WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
   WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
   WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
   WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
   WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
   WRITE(IOUNLOG,*)
   WRITE(IOUNLOG,*)

   !... EXPLICIT VERTICAL DIFFUSION
   CALL RANDOM_FLD(GTEM)
   CALL RANDOM_FLD(GAS)
   CALL RANDOM_FLD(GAT)
   CALL RANDOM_FLD(GTEM_AD)
   GAT0=GAT
   CALL TRA_ZDF_EXP(KN_ZDFEXP,GTEM,GAT,AVBT)
   ZT2 = MDOT_PRODUCT(GAT,GAS)
   CALL TRA_ZDF_EXP_AD(KN_ZDFEXP,GTEM_AD,GAS,AVBT)
   ZT1 = MDOT_PRODUCT(GTEM_AD,GTEM) + MDOT_PRODUCT(GAT0,GAS)
   WRITE(IOUNLOG,*) ' ADJOINT TEST FOR TRA_ZDF_EXP'
   WRITE(IOUNLOG,*) ' GTEM EXTRMS = ',MINVAL(GTEM),MAXVAL(GTEM),&
   &SUM(GTEM)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAS EXTRMS = ',MINVAL(GAS),MAXVAL(GAS),&
   &SUM(GAS)/SUM(GRD%MSK) 
   WRITE(IOUNLOG,*) ' GTEM_AD EXTRMS = ',MINVAL(GTEM_AD),MAXVAL(GTEM_AD),&
   &SUM(GTEM_AD)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' GAT EXTRMS = ',MINVAL(GAT),MAXVAL(GAT),&
   &SUM(GAT)/SUM(GRD%MSK)
   WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
   WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
   WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
   WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
   WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
   WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
   WRITE(IOUNLOG,*)
   WRITE(IOUNLOG,*)

   !... SURFACE FLUXES
   CALL RANDOM_FLD(GTEM2)
   CALL RANDOM_FLD(GSAL2)
   CALL RANDOM_FLD(GAT2B)
   CALL RANDOM_FLD(GAS2B)
   GAT2 = 0._R8
   GAS2 = 0._R8
   GTEM2_AD = 0._R8
   GSAL2_AD = 0._R8
   GAT2B0 = GAT2B
   GAS2B0 = GAS2B
   WRITE(IOUNLOG,*) 'GTEM2',MINVAL(GTEM2),MAXVAL(GTEM2)
   WRITE(IOUNLOG,*) 'GSAL2',MINVAL(GSAL2),MAXVAL(GSAL2)
   CALL BLK_OCE_CORE(JPI,JPJ,GTEM2,GSAL2,GAT2,GAS2)
   CALL BLK_OCE_CORE_AD(JPI,JPJ,GTEM2_AD,GSAL2_AD,GAT2B,GAS2B)
   WRITE(IOUNLOG,*) 'GTEM2',MINVAL(GTEM2),MAXVAL(GTEM2)
   WRITE(IOUNLOG,*) 'GSAL2',MINVAL(GSAL2),MAXVAL(GSAL2)
   ALLOCATE( ZV1(GRD%IM*GRD%JM*2) )
   ALLOCATE( ZV2(GRD%IM*GRD%JM*2) )
   ALLOCATE( ZV3(GRD%IM*GRD%JM*2) )
   ALLOCATE( ZV4(GRD%IM*GRD%JM*2) )
   KK=0
   DO JJ=1,GRD%JM
     DO JI=1,GRD%IM
       KK=KK+1
       ZV1(KK) = GTEM2(JI,JJ)
       ZV2(KK) = GTEM2_AD(JI,JJ)
       ZV3(KK) = GAT2(JI,JJ)
       ZV4(KK) = GAT2B0(JI,JJ)
     ENDDO
   ENDDO
   DO JJ=1,GRD%JM
     DO JI=1,GRD%IM
       KK=KK+1
       ZV1(KK) = GSAL2(JI,JJ)
       ZV2(KK) = GSAL2_AD(JI,JJ)
       ZV3(KK) = GAS2(JI,JJ)
       ZV4(KK) = GAS2B0(JI,JJ)
     ENDDO
   ENDDO
   ZT1 = DOT_PRODUCT(ZV1,ZV2)
   ZT2 = DOT_PRODUCT(ZV3,ZV4)
   WRITE(IOUNLOG,*) ' ADJOINT TEST FOR SURFACE FLUXES'
   WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
   WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
   WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
   WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
   WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
   WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
   WRITE(IOUNLOG,*)
   WRITE(IOUNLOG,*)
   DEALLOCATE(ZV1,ZV2,ZV3,ZV4)

   CALL RANDOM_FLD(GRD%TEM)
   CALL RANDOM_FLD(GRD%SAL)
   GAT0=GRD%TEM
   GAS0=GRD%SAL
   ALLOCATE( ZV1(GRD%IM*GRD%JM*GRD%KM*2) )
   ALLOCATE( ZV2(GRD%IM*GRD%JM*GRD%KM*2) )
   ALLOCATE( ZV3(GRD%IM*GRD%JM*GRD%KM*2) )
   ALLOCATE( ZV4(GRD%IM*GRD%JM*GRD%KM*2) )
   LLSAVE=LL_WEAKLY
   LL_WEAKLY=.FALSE.
   CALL TLMOD(2)
   CALL RANDOM_FLD(GRD%TEM_AD)
   CALL RANDOM_FLD(GRD%SAL_AD)
   GAT=GRD%TEM_AD
   GAS=GRD%SAL_AD
   CALL ADMOD(2)
   KK=0
   DO JK=1,GRD%KM
     DO JJ=1,GRD%JM
       DO JI=1,GRD%IM
         KK=KK+1
         ZV1(KK) = GAT0(JI,JJ,JK)
         ZV2(KK) = GRD%TEM_AD(JI,JJ,JK)
         ZV3(KK) = GRD%TEM(JI,JJ,JK)
         ZV4(KK) = GAT(JI,JJ,JK)
       ENDDO
     ENDDO
   ENDDO
   DO JK=1,GRD%KM
     DO JJ=1,GRD%JM
       DO JI=1,GRD%IM
         KK=KK+1
         ZV1(KK) = GAS0(JI,JJ,JK)
         ZV2(KK) = GRD%SAL_AD(JI,JJ,JK)
         ZV3(KK) = GRD%SAL(JI,JJ,JK)
         ZV4(KK) = GAS(JI,JJ,JK)
       ENDDO
     ENDDO
   ENDDO
   ZT1 = DOT_PRODUCT(ZV1,ZV2)
   ZT2 = DOT_PRODUCT(ZV3,ZV4)
   WRITE(IOUNLOG,*) ' ADJOINT TEST FOR TLMOD'
   WRITE(IOUNLOG,*) ' FIRST  TERM = ',ZT1
   WRITE(IOUNLOG,*) ' SECOND TERM = ',ZT2
   WRITE(IOUNLOG,*) ' DIFFERENCE  = ',ZT1-ZT2
   WRITE(IOUNLOG,*) ' MACHINE EPS = ',EPSILON(ZT1)
   WRITE(IOUNLOG,*) ' REL. ERROR  = ',100._R8*ABS((ZT1-ZT2)/ZT1)
   WRITE(IOUNLOG,*) ' RESULT      = ',TEST_RES(ABS(ZT1-ZT2),EPSILON(ZT1))
   WRITE(IOUNLOG,*)
   WRITE(IOUNLOG,*)
   
   LL_WEAKLY = LLSAVE
   GRD%TEM=0._R8
   GRD%SAL=0._R8
   GRD%TEM_AD=0._R8
   GRD%TEM_AD=0._R8
   DEALLOCATE(ZV1,ZV2,ZV3,ZV4)

 CALL MYFRTPROF_WALL('TLMOD_TEST: TEST ADJOINT',1)

END SUBROUTINE TLADMOD_TEST

END MODULE TLAD_TEST
