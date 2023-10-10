MODULE TRAZDF_IMP

   USE TLAD_VARS
   USE MYFRTPROF

   IMPLICIT NONE

CONTAINS

   SUBROUTINE TRA_ZDF_IMP(RDT2, PTB_TL, PTA_TL, AVT)
   IMPLICIT NONE
      REAL(WP), INTENT(IN   ) ::   RDT2
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(IN   ) ::   PTB_TL   ! BEFORE AND NOW TRACER FIELDS
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(IN   ) ::   AVT  
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(INOUT) ::   PTA_TL   ! TRACER TREND
      !! * LOCAL DECLARATIONS
      INTEGER  ::   JI, JJ, JK               ! DUMMY LOOP INDICES
      REAL(WP) ::   ZRHSTL,          & ! TEMPORARY SCALARS
         ZE3TB, ZE3TN, ZE3TA, ZVSFVVL        ! VARIABLE VERTICAL SCALE FACTORS
      REAL(WP), POINTER, DIMENSION(:,:,:) ::   &
         ZWI, ZWT, ZWD, ZWS                     ! WORKSPACE ARRAYS
      !!---------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('TRA_ZDF_IMP: TL IMPLICIT VERT. MIXING',0)

      ALLOCATE( ZWI(JPI,JPJ,JPK) )
      ALLOCATE( ZWT(JPI,JPJ,JPK) )
      ALLOCATE( ZWD(JPI,JPJ,JPK) )
      ALLOCATE( ZWS(JPI,JPJ,JPK) )
      ZWI = 0._WP
      ZWT = 0._WP
      ZWD = 0._WP
      ZWS = 0._WP
      !
         ! --------------------
         ! BUILD MATRIX IF TEMPERATURE OR SALINITY (ONLY IN DOUBLE DIFFUSION CASE) OR FIRST PASSIVE TRACER
         !
            ZWT(:,:,2:JPK) = AVT  (:,:,2:JPK)
            ZWT(:,:,1) = 0._WP
            !
            ! II.0 MATRIX CONSTRUCTION
            ! ------------------------
            ! DIAGONAL, INFERIOR, SUPERIOR  (INCLUDING THE BOTTOM BOUNDARY CONDITION VIA AVT MASKED)
!$omp parallel  default(shared) private(jj,jk,ji,ze3ta,ze3tb,ze3tn,zrhstl)
!$omp do schedule(dynamic,1)
               DO JJ = 2, JPJM1
            DO JK = 1, JPKM1
                  DO JI = 2, JPIM1   ! VECTOR OPT.
                     ZE3TA = 1._WP                                ! AFTER SCALE FACTOR AT T-POINT
                     ZE3TN = E3T(JK)                      ! NOW   SCALE FACTOR AT T-POINT
                     ZWI(JI,JJ,JK) = - RDT2 * ZWT(JI,JJ,JK  ) / ( ZE3TN * E3W(JK  ) )
                     ZWS(JI,JJ,JK) = - RDT2 * ZWT(JI,JJ,JK+1) / ( ZE3TN * E3W(JK+1) )
                     ZWD(JI,JJ,JK) = ZE3TA - ZWI(JI,JJ,JK) - ZWS(JI,JJ,JK)
                  END DO
               END DO
!            END DO

            ! FIRST RECURRENCE:   TK = DK - IK SK-1 / TK-1   (INCREASING K)
!            DO JJ = 2, JPJM1
               DO JI = 2, JPIM1
                  ZWT(JI,JJ,1) = ZWD(JI,JJ,1)
               END DO
!            END DO
            DO JK = 2, JPKM1
!               DO JJ = 2, JPJM1
                  DO JI = 2, JPIM1
                     ZWT(JI,JJ,JK) = ZWD(JI,JJ,JK) - ZWI(JI,JJ,JK) * ZWS(JI,JJ,JK-1)  /ZWT(JI,JJ,JK-1)
                  END DO
!               END DO
            END DO

         ! SECOND RECURRENCE:    ZK = YK - IK / TK-1  ZK-1
!         DO JJ = 2, JPJM1
            DO JI = 2, JPIM1
               ZE3TB = 1._WP
               ZE3TN = 1._WP
               PTA_TL(JI,JJ,1) = ZE3TB * PTB_TL(JI,JJ,1) + RDT2 * ZE3TN * PTA_TL(JI,JJ,1)
            END DO
!         END DO

         DO JK = 2, JPKM1
!            DO JJ = 2, JPJM1
               DO JI = 2, JPIM1
                  ZE3TB = 1._WP
                  ZE3TN = 1._WP
                  ZRHSTL = ZE3TB * PTB_TL(JI,JJ,JK) + RDT2 * ZE3TN * PTA_TL(JI,JJ,JK)   ! ZRHS=RIGHT HAND SIDE
                  PTA_TL(JI,JJ,JK) = ZRHSTL - ZWI(JI,JJ,JK) / ZWT(JI,JJ,JK-1) * PTA_TL(JI,JJ,JK-1)
               END DO
!            END DO
         END DO

         ! THIRD RECURRENCE: XK = (ZK - SK XK+1 ) / TK
         ! SAVE THE MASKED TEMPERATURE AFTER IN TA
         ! (C A U T I O N:  TEMPERATURE NOT ITS TREND, LEAP-FROG SCHEME DONE IT WILL NOT BE DONE IN TRANXT)
!         DO JJ = 2, JPJM1
            DO JI = 2, JPIM1
               PTA_TL(JI,JJ,JPKM1) = PTA_TL(JI,JJ,JPKM1) / ZWT(JI,JJ,JPKM1) * TMASK(JI,JJ,JPKM1)
!            END DO
         END DO
         DO JK = JPK-2, 1, -1
!            DO JJ = 2, JPJM1
               DO JI = 2, JPIM1
                  PTA_TL(JI,JJ,JK) = ( PTA_TL(JI,JJ,JK) - ZWS(JI,JJ,JK) * PTA_TL(JI,JJ,JK+1) )   &
                                  &   / ZWT(JI,JJ,JK) * TMASK(JI,JJ,JK)
               END DO
            END DO
         END DO
! END of openMP loop on JJ
!$omp end do
!$omp end parallel
      !
      DEALLOCATE( ZWI, ZWT, ZWD, ZWS )
      CALL MYFRTPROF_WALL('TRA_ZDF_IMP: TL IMPLICIT VERT. MIXING',1)
      !
   END SUBROUTINE TRA_ZDF_IMP

!-----------------------------------------------------------------------

   SUBROUTINE TRA_ZDF_IMP_AD(RDT2,PTB_AD, PTA_AD, AVT)
   IMPLICIT NONE
      REAL(WP), INTENT(IN   ) ::   RDT2
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(IN   ) ::   AVT  
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(INOUT) ::   PTB_AD      ! BEFORE AND NOW TRACER FIELDS
      REAL(WP), DIMENSION(JPI,JPJ,JPK), INTENT(INOUT) ::   PTA_AD      ! TRACER TREND
      !! * LOCAL DECLARATIONS
      INTEGER  ::   JI, JJ, JK               ! DUMMY LOOP INDICES
      REAL(WP) ::   ZTEMPB
      REAL(WP) ::   ZRHSAD,          & ! TEMPORARY SCALARS
         ZE3TB, ZE3TN, ZE3TA, ZVSFVVL        ! VARIABLE VERTICAL SCALE FACTORS
      REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) ::   &
         ZWI, ZWT, ZWS, ZWD                     ! WORKSPACE ARRAYS
      !!---------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('TRA_ZDF_IMP_AD: AD IMPLICIT VERT. MIXING',0)
      ALLOCATE ( ZWI(JPI,JPJ,JPK) )
      ALLOCATE ( ZWT(JPI,JPJ,JPK) )
      ALLOCATE ( ZWS(JPI,JPJ,JPK) )
      ALLOCATE ( ZWD(JPI,JPJ,JPK) )
      ZWI = 0._WP
      ZWT = 0._WP
      ZWD = 0._WP
      ZWS = 0._WP
      !
            ZWT(:,:,2:JPK) = AVT  (:,:,2:JPK)
            ZWT(:,:,1) = 0._WP

            ! DIAGONAL, INFERIOR, SUPERIOR  (INCLUDING THE BOTTOM BOUNDARY CONDITION VIA AVT MASKED)
!$omp parallel  default(shared) private(jj,jk,ji,ze3ta,ze3tn,ze3tb,zrhsad,ztempb)
!$omp do schedule(dynamic,1)
               DO JJ = 2, JPJM1
            DO JK = 1, JPKM1
                  DO JI = 2, JPIM1   ! VECTOR OPT.
                     ZE3TA = 1._WP                                ! AFTER SCALE FACTOR AT T-POINT
                     ZE3TN = E3T(JK)                      ! NOW   SCALE FACTOR AT T-POINT
                     ZWI(JI,JJ,JK) = - RDT2 * ZWT(JI,JJ,JK  ) / ( ZE3TN * E3W(JK  ) )
                     ZWS(JI,JJ,JK) = - RDT2 * ZWT(JI,JJ,JK+1) / ( ZE3TN * E3W(JK+1) )
                     ZWD(JI,JJ,JK) = ZE3TA - ZWI(JI,JJ,JK) - ZWS(JI,JJ,JK)
                  END DO
!               END DO
            END DO

            !! MATRIX INVERSION FROM THE FIRST LEVEL
            !!----------------------------------------------------------------------
!            DO JJ = 2, JPJM1
               DO JI = 2, JPIM1
                  ZWT(JI,JJ,1) = ZWD(JI,JJ,1)
               END DO
!            END DO
            DO JK = 2, JPKM1
!               DO JJ = 2, JPJM1
                  DO JI = 2, JPIM1
                     ZWT(JI,JJ,JK) = ZWD(JI,JJ,JK) - ZWI(JI,JJ,JK) * ZWS(JI,JJ,JK-1)  /ZWT(JI,JJ,JK-1)
                  END DO
!               END DO
            END DO

         ! THIRD RECURRENCE: XK = (ZK - SK XK+1 ) / TK
         DO JK = 1, JPK-2
!            DO JJ = JPJM1,2,-1
               DO JI = JPIM1,2,-1
                  PTA_AD(JI,JJ,JK+1) = PTA_AD(JI,JJ,JK+1) - ZWS(JI,JJ,JK) * PTA_AD(JI,JJ,JK)   &
                                    &   / ZWT(JI,JJ,JK) * TMASK(JI,JJ,JK)
                  PTA_AD(JI,JJ,JK)   = PTA_AD(JI,JJ,JK) / ZWT(JI,JJ,JK) * TMASK(JI,JJ,JK)
               END DO
!            END DO
         END DO
!         DO JJ = JPJM1,2,-1
            DO JI = JPIM1,2,-1
               PTA_AD(JI,JJ,JPKM1) = TMASK(JI,JJ,JPKM1)*PTA_AD(JI,JJ,JPKM1) / ZWT(JI,JJ,JPKM1)
            END DO
!         END DO
         ! SECOND RECURRENCE:    ZK = YK - IK / TK-1  ZK-1
!.....         PTB_AD(:,JJ,:) = 0._WP
         DO JK = JPKM1, 2, -1
!            DO JJ = JPJM1,2,-1
               DO JI = JPIM1,2,-1
                  ZE3TB = 1._WP
                  ZE3TN = 1._WP
                  ZRHSAD = PTA_AD(JI,JJ,JK)
                  PTA_AD(JI,JJ,JK-1) = PTA_AD(JI,JJ,JK-1)-ZWI(JI,JJ,JK)*PTA_AD(JI,JJ,JK)/ZWT(JI,JJ,JK-1)
                  PTA_AD(JI,JJ,JK) = 0._WP
                  PTB_AD(JI,JJ,JK) = PTB_AD(JI,JJ,JK) + ZE3TB * ZRHSAD
                  PTA_AD(JI,JJ,JK) = PTA_AD(JI,JJ,JK) + RDT2 * ZE3TN * ZRHSAD
               END DO
!            END DO
         END DO
!         DO JJ = JPJM1,2,-1
            DO JI = JPIM1,2,-1
               ZE3TB = 1._WP
               ZE3TN = 1._WP
               PTB_AD(JI,JJ,1) = PTB_AD(JI,JJ,1) + ZE3TB * PTA_AD(JI,JJ,1)
               PTA_AD(JI,JJ,1) = PTA_AD(JI,JJ,1) * RDT2 * ZE3TN
            END DO
         END DO
!$omp end do
!$omp end parallel
      !
      DEALLOCATE( ZWI, ZWT, ZWS, ZWD )
      CALL MYFRTPROF_WALL('TRA_ZDF_IMP_AD: AD IMPLICIT VERT. MIXING',1)
      !
   END SUBROUTINE TRA_ZDF_IMP_AD
END MODULE TRAZDF_IMP
