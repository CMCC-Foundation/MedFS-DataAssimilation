MODULE traadv_cen2h

   USE TLAD_VARS
   USE MYFRTPROF
   IMPLICIT NONE

   PRIVATE
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwx, zwy
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwxad, zwyad
   LOGICAL, PRIVATE, SAVE :: LLINIT_TL = .FALSE.
   LOGICAL, PRIVATE, SAVE :: LLINIT_AD = .FALSE.

   PUBLIC tra_adv_cen2h, tra_adv_cen2h_ad

  INTERFACE tra_adv_cen2h
         MODULE PROCEDURE  tra_adv_cen2_3d
  END INTERFACE

  INTERFACE tra_adv_cen2h_ad
         MODULE PROCEDURE  tra_adv_cen2_ad_3d
  END INTERFACE

CONTAINS

   SUBROUTINE tra_adv_cen2_3d(jpi,jpj,jpk,pun,pvn,ptn,pta)

      IMPLICIT NONE

      !
      INTEGER, INTENT( in ) ::   jpi,jpj,jpk
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: pun,pvn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: ptn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(INOUT) :: pta

      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ierr             ! local integer
      REAL(wp) ::   zbtr, ztra                            ! local scalars
      REAL(wp) ::   zcenut, zcenvt, zcent
      !!----------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('TRA_ADV_CEN2H_3D: TL HORIZ ADVECTION',0)

      IF( .NOT. LLINIT_TL ) THEN
        ALLOCATE( zwx(jpi,jpj,jpk), &
                zwy(jpi,jpj,jpk) )
        zwx = 0._wp
        zwy = 0._wp
        LLINIT_TL = .TRUE.
      ENDIF
      !                                          ! set time step
      !
      ! I. Horizontal advection
      !    ====================
      !
!$omp parallel  default(shared) private(jk,jj,ji,zcenut,zcenvt)
!$omp do schedule(dynamic,1)
      DO jk = 1, jpkm1
         ! Second order centered tracer flux at u- and v-points
         DO jj = 1, jpjm1
            !
            DO ji = 1, jpim1   ! vector opt.
               ! centered scheme
               zcenut = pun(ji,jj,jk) * ( ptn(ji,jj,jk) + ptn(ji+1,jj  ,jk) ) * e2u(ji,jj) 
               zcenvt = pvn(ji,jj,jk) * ( ptn(ji,jj,jk) + ptn(ji  ,jj+1,jk) ) * e1v(ji,jj)
               ! mixed centered / upstream scheme
               zwx(ji,jj,jk) = 0.5 * zcenut
               zwy(ji,jj,jk) = 0.5 * zcenvt 
            END DO
         END DO
      END DO
!$omp end do
!$omp end parallel

         ! II. Divergence of advective fluxes
         ! ----------------------------------
!$omp parallel  default(shared) private(jk,jj,ji,zbtr,ztra)
!$omp do schedule(dynamic,1)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) )
               ! advective trends
               ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
               &                + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   )
               ! advective trends added to the general tracer trends
               pta(ji,jj,jk) = pta(ji,jj,jk) + ztra
            END DO
         END DO
      END DO
!$omp end do
!$omp end parallel

      CALL MYFRTPROF_WALL('TRA_ADV_CEN2H_3D: TL HORIZ ADVECTION',1)

   END SUBROUTINE tra_adv_cen2_3d

   SUBROUTINE tra_adv_cen2_ad_3d(jpi,jpj,jpk,pun,pvn,ptn_ad,pta_ad)

      IMPLICIT NONE

      !
      INTEGER, INTENT( in ) ::   jpi,jpj,jpk
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: pun,pvn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(INOUT) :: ptn_ad
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: pta_ad

      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   ierr             ! local integer
      REAL(wp) ::   zbtr, ztra                            ! local scalars
      REAL(wp) ::   zfui, zfvj, zhw
      !!----------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('TRA_ADV_CEN2H_AD_3D: AD HORIZ ADVECTION',0)

      !
      IF( .NOT. LLINIT_AD ) THEN
        ALLOCATE( zwxad(jpi,jpj,jpk), &
                zwyad(jpi,jpj,jpk) )
        LLINIT_AD = .TRUE.
      ENDIF
      zwxad = 0._R8
      zwyad = 0._R8
      !                                          ! set time step
      ! II. Divergence of advective fluxes
      ! ----------------------------------
!$omp parallel  default(shared) private(jk,jj,ji,zbtr)
!$omp do schedule(dynamic,1)
      DO jk = jpkm1, 1, -1
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) )
               zwxad(ji,jj,jk) = zwxad(ji,jj,jk) - zbtr * pta_ad(ji,jj,jk)
               zwyad(ji,jj,jk) = zwyad(ji,jj,jk) - zbtr * pta_ad(ji,jj,jk)
               zwxad(ji-1,jj,jk) = zwxad(ji-1,jj,jk) + zbtr * pta_ad(ji,jj,jk)
               zwyad(ji,jj-1,jk) = zwyad(ji,jj-1,jk) + zbtr * pta_ad(ji,jj,jk)
            END DO
         END DO
      END DO
!$omp end do
!$omp end parallel

!$omp parallel  default(shared) private(jk,jj,ji,zfui,zfvj)
!$omp do schedule(dynamic,1)
      DO jk = jpkm1, 1, -1
         ! Second order centered tracer flux at u- and v-points
         DO jj = jpjm1, 1, -1
            DO ji = jpim1, 1, -1   ! vector opt.
               ! volume fluxes * 1/2
               zfvj = 0.5 * pvn(ji,jj,jk)*zwyad(ji,jj,jk)*e1v(ji,jj)
               zwyad(ji,jj,jk) = 0._wp
               zfui = 0.5 * pun(ji,jj,jk)*zwxad(ji,jj,jk)*e2u(ji,jj)
               zwxad(ji,jj,jk) = 0._wp
               ! centered scheme
               ptn_ad(ji  ,jj  ,jk) = ptn_ad(ji  ,jj  ,jk) + zfvj
               ptn_ad(ji  ,jj+1,jk) = ptn_ad(ji  ,jj+1,jk) + zfvj
               ptn_ad(ji  ,jj  ,jk) = ptn_ad(ji  ,jj  ,jk) + zfui
               ptn_ad(ji+1,jj  ,jk) = ptn_ad(ji+1,jj  ,jk) + zfui
            END DO
         END DO
      END DO
!$omp end do
!$omp end parallel

      CALL MYFRTPROF_WALL('TRA_ADV_CEN2H_AD_3D: AD HORIZ ADVECTION',1)

   END SUBROUTINE tra_adv_cen2_ad_3d

END MODULE traadv_cen2h
