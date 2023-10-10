MODULE TRAZDF_EXP

   USE TLAD_VARS
   USE IOUNITS
   USE MYFRTPROF

   IMPLICIT NONE

CONTAINS

   SUBROUTINE TRA_ZDF_EXP(KN_ZDFEXP,PTB,PTA,AVT)
    
      INTEGER                              , INTENT(IN   ) ::   KN_ZDFEXP   ! NUMBER OF SUB-TIME STEP
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(IN   ) ::   PTB         ! BEFORE AND NOW TRACER FIELDS
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(INOUT) ::   PTA         ! TRACER TREND
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(IN   ) ::   AVT 
      !!
      INTEGER  ::   JI, JJ, JK, JL            ! DUMMY LOOP INDICES
      REAL(WP) ::   ZLAVMR, ZAVE3R, ZE3TR     ! TEMPORARY SCALARS
      REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) ::   ZWXTL, ZWYTL   ! 3D WORKSPACE
      !!---------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('TRA_ZDF_EXP: TL EXPLICIT VERT. MIXING',0)
      ALLOCATE( ZWXTL(JPI, JPJ, JPK) )
      ALLOCATE( ZWYTL(JPI, JPJ, JPK) )
      !
      ! INITIALIZATIONS
      ! ---------------
      ZLAVMR = 1. / FLOAT( KN_ZDFEXP )                           ! LOCAL CONSTANT
      !
         ZWYTL(:,:, 1 ) = 0.0_WP                ! SURFACE BOUNDARY CONDITIONS: NO FLUX
         ZWYTL(:,:,JPK) = 0.0_WP                ! BOTTOM  BOUNDARY CONDITIONS: NO FLUX
         !
         ZWXTL(:,:,:)   = PTB(:,:,:)      ! ZWX AND ZWZ ARRAYS SET TO BEFORE TRACER VALUES

         ! SPLIT-EXPLICIT LOOP  (AFTER TRACER DUE TO THE VERTICAL DIFFUSION ALONE)
         ! -------------------
         !
         IF(LL_TLDEB) WRITE(IOUNLOG,*) ' ** TS',0,ZWXTL(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1)),&
         & ZWYTL(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1)),AVT(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1))
         DO JL = 1, KN_ZDFEXP
            !                     ! FIRST VERTICAL DERIVATIVE
            DO JK = 2, JPK
               DO JJ = 2, JPJM1
                  DO JI = 2, JPIM1   ! VECTOR OPT.
                     ZAVE3R = 1._WP / E3W(JK)
                     ZWYTL(JI,JJ,JK) = AVT(JI,JJ,JK) * (ZWXTL(JI,JJ,JK-1)-ZWXTL(JI,JJ,JK)) &
                     & * ZAVE3R
                  END DO
               END DO
            END DO
            !
            DO JK = 1, JPKM1   ! SECOND VERTICAL DERIVATIVE ==> TRACER AT KT+L*2*RDT/N_ZDFEXP
               DO JJ = 2, JPJM1
                  DO JI = 2, JPIM1   ! VECTOR OPT.
                     ZE3TR = ZLAVMR / E3T(JK)
                     ZWXTL(JI,JJ,JK) = ZWXTL(JI,JJ,JK) + &
                     & RDT*(ZWYTL(JI,JJ,JK)-ZWYTL(JI,JJ,JK+1) ) * ZE3TR
                  END DO
               END DO
            END DO
            IF(LL_TLDEB) WRITE(IOUNLOG,*) ' ** TS',JL,ZWXTL(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1)),&
            & ZWYTL(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1)),AVT(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1))

            !
         END DO

         ! AFTER TRACER DUE TO ALL TRENDS
         ! ------------------------------
         IF(LL_TLDEB) WRITE(IOUNLOG,*) ' ** TS',-1,ZWXTL(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1)),PTA(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1))
         DO JK = 1, JPKM1
            DO JJ = 2, JPJM1
               DO JI = 2, JPIM1   ! VECTOR OPT.
                  PTA(JI,JJ,JK) = ( ZWXTL(JI,JJ,JK) + RDT * PTA(JI,JJ,JK) ) * TMASK(JI,JJ,JK)
               END DO
            END DO
         END DO
         IF(LL_TLDEB) WRITE(IOUNLOG,*) ' ** TS',-2,ZWXTL(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1)),PTA(IDB_COO(1),IDB_COO(2),(IDB_COO(3)-1):(IDB_COO(3)+1))
      !
      DEALLOCATE( ZWXTL, ZWYTL)
      CALL MYFRTPROF_WALL('TRA_ZDF_EXP: TL EXPLICIT VERT. MIXING',1)
      !
   END SUBROUTINE TRA_ZDF_EXP

   SUBROUTINE TRA_ZDF_EXP_AD(KN_ZDFEXP,PTB_AD,PTA_AD,AVT)
      INTEGER                              , INTENT(IN   ) ::   KN_ZDFEXP   ! NUMBER OF SUB-TIME STEP
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(INOUT) ::   PTB_AD      ! BEFORE AND NOW TRACER FIELDS
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(INOUT) ::   PTA_AD      ! TRACER TREND
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(IN   ) ::   AVT 
      !!
      INTEGER  ::   JI, JJ, JK, JL            ! DUMMY LOOP INDICES
      REAL(WP) ::   ZLAVMR, ZAVE3R, ZE3TR     ! TEMPORARY SCALARS
      REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) ::   ZWXAD, ZWYAD                 ! 3D WORKSPACE
      !!---------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('TRA_ZDF_EXP_AD: AD EXPLICIT VERT. MIXING',0)
      ALLOCATE( ZWXAD(JPI, JPJ, JPK) )
      ALLOCATE( ZWYAD(JPI, JPJ, JPK) )
      !
      ! INITIALIZATIONS
      ! ---------------
      ZLAVMR = 1. / FLOAT( KN_ZDFEXP )                           ! LOCAL CONSTANT
         !
         ZWXAD(:,:,:) = 0.0_WP
         ZWYAD(:,:,:) = 0.0_WP
         ! AFTER TRACER DUE TO ALL TRENDS
         ! ------------------------------
         DO JK = 1, JPKM1
            DO JJ = 2, JPJM1
               DO JI = 2, JPIM1   ! VECTOR OPT.
                  ZWXAD(JI,JJ,JK) = ZWXAD(JI,JJ,JK) + PTA_AD(JI,JJ,JK) * TMASK(JI,JJ,JK)
                  PTA_AD(JI,JJ,JK) = RDT * PTA_AD(JI,JJ,JK) * TMASK(JI,JJ,JK)
               END DO
            END DO
         END DO
         !

         ! SPLIT-EXPLICIT LOOP  (AFTER TRACER DUE TO THE VERTICAL DIFFUSION ALONE)
         ! -------------------
         !
         DO JL = 1, KN_ZDFEXP
            DO JK =  JPKM1, 1, -1      ! SECOND VERTICAL DERIVATIVE   ==> TRACER AT KT+L*2*RDT/N_ZDFEXP
               DO JJ = 2, JPJM1
                  DO JI = 2, JPIM1   ! VECTOR OPT.
                     ZE3TR = ZLAVMR / E3T(JK)
                     ZWYAD(JI,JJ,JK  ) = ZWYAD(JI,JJ,JK  ) + RDT * ZWXAD(JI,JJ,JK) * ZE3TR
                     ZWYAD(JI,JJ,JK+1) = ZWYAD(JI,JJ,JK+1) - RDT * ZWXAD(JI,JJ,JK) * ZE3TR
                  END DO
               END DO
            END DO
            !                     ! FIRST VERTICAL DERIVATIVE
            DO JK = JPK, 2, -1
               DO JJ = 2, JPJM1
                  DO JI = 2, JPIM1   ! VECTOR OPT.
                     ZAVE3R = 1._WP / E3W(JK)
                     ZWXAD(JI,JJ,JK-1) = ZWXAD(JI,JJ,JK-1) + AVT(JI,JJ,JK) * &
                     & ZWYAD(JI,JJ,JK) * ZAVE3R
                     ZWXAD(JI,JJ,JK  ) = ZWXAD(JI,JJ,JK  ) - AVT(JI,JJ,JK) * &
                     & ZWYAD(JI,JJ,JK) * ZAVE3R
                     ZWYAD(JI,JJ,JK  ) = 0.0_WP
                  END DO
               END DO
            END DO
            !
            !
         END DO
         !
         PTB_AD(:,:,:) = PTB_AD(:,:,:) + ZWXAD(:,:,:)
      DEALLOCATE( ZWXAD, ZWYAD )
      CALL MYFRTPROF_WALL('TRA_ZDF_EXP_AD: AD EXPLICIT VERT. MIXING',1)
      !
   END SUBROUTINE TRA_ZDF_EXP_AD

END MODULE TRAZDF_EXP
