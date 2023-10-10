MODULE traadv_cen2v

   USE TLAD_VARS
   USE MYFRTPROF
   IMPLICIT NONE

   PRIVATE 
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwz
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwzad
   LOGICAL, PRIVATE, SAVE :: LLINIT_TL = .FALSE.
   LOGICAL, PRIVATE, SAVE :: LLINIT_AD = .FALSE.

   PUBLIC tra_adv_cen2v, tra_adv_cen2v_ad

  INTERFACE tra_adv_cen2v
         MODULE PROCEDURE  tra_adv_cen2_3d
  END INTERFACE

  INTERFACE tra_adv_cen2v_ad
         MODULE PROCEDURE  tra_adv_cen2_ad_3d
  END INTERFACE

CONTAINS

   SUBROUTINE tra_adv_cen2_3d(jpi,jpj,jpk,pwn,ptn,pta)

      IMPLICIT NONE

      !
      INTEGER, INTENT( in ) ::   jpi,jpj,jpk
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: pwn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: ptn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(INOUT) :: pta

      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ierr             ! local integer
      REAL(wp) ::   zbtr, ztra                            ! local scalars
      REAL(wp) ::   zcenut, zcenvt, zcent
      !!----------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('TRA_ADV_CEN2V_3D: TL VERT ADVECTION',0)

      IF( .NOT. LLINIT_TL ) THEN
        ALLOCATE( zwz(jpi,jpj,jpk)  )
        zwz = 0._wp
        LLINIT_TL = .TRUE.
      ENDIF
      !                                          ! set time step
      zwz(:,:,jpk) = 0.e0             ! Bottom  value : flux set to zero
      ! Surface value : 
      zwz(:,:, 1 ) = pwn(:,:,1) * ptn(:,:,1)   ! linear free surface 
      !
!$omp parallel  default(shared) private(jk,jj,ji,zcent)
!$omp do schedule(dynamic,1)
      DO jk = 2, jpk              ! Second order centered tracer flux at w-point
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ! centered scheme
               zcent = pwn(ji,jj,jk) * ( ptn(ji,jj,jk) + ptn(ji,jj,jk-1) )
               ! mixed centered / upstream scheme
               zwz(ji,jj,jk) = 0.5 * zcent
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
               zbtr = 1./e3t(jk)
               ! advective trends
               ztra = - zbtr * ( zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1)  )
               ! advective trends added to the general tracer trends
               pta(ji,jj,jk) = pta(ji,jj,jk) + ztra
            END DO
         END DO
      END DO
!$omp end do
!$omp end parallel

      CALL MYFRTPROF_WALL('TRA_ADV_CEN2V_3D: TL VERT ADVECTION',1)

   END SUBROUTINE tra_adv_cen2_3d

   SUBROUTINE tra_adv_cen2_ad_3d(jpi,jpj,jpk,pwn,ptn_ad,pta_ad)

      IMPLICIT NONE

      !
      INTEGER, INTENT( in ) ::   jpi,jpj,jpk
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: pwn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(INOUT) :: ptn_ad
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: pta_ad

      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   ierr             ! local integer
      REAL(wp) ::   zbtr, ztra                            ! local scalars
      REAL(wp) ::   zfui, zfvj, zhw
      !!----------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('TRA_ADV_CEN2V_AD_3D: AD VERT ADVECTION',0)

      !
      IF( .NOT. LLINIT_AD ) THEN
        ALLOCATE( zwzad(jpi,jpj,jpk)  )
        LLINIT_AD = .TRUE.
      ENDIF
      zwzad = 0._R8
      !                                          ! set time step
!$omp parallel  default(shared) private(jk,jj,ji,zbtr)
!$omp do schedule(dynamic,1)
         DO jj = jpjm1, 2, -1
      DO jk = jpkm1, 1, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               zbtr = 1. / ( e3t(jk) )
               zwzad(ji,jj,jk) = zwzad(ji,jj,jk) - zbtr * pta_ad(ji,jj,jk)
               zwzad(ji,jj,jk+1) = zwzad(ji,jj,jk+1) + zbtr * pta_ad(ji,jj,jk)
            END DO
         END DO
      END DO
!$omp end do
!$omp end parallel

!$omp parallel  default(shared) private(jj,jk,ji,zhw)
!$omp do schedule(dynamic,1)
         DO jj = jpjm1, 2, -1
      DO jk = jpk, 2, -1      ! Second order centered tracer flux at w-point
            DO ji = jpim1, 2, -1   ! vector opt.
               zhw   = 0.5 * pwn(ji,jj,jk)*zwzad(ji,jj,jk)
               zwzad(ji,jj,jk) = 0._R8
               ! centered scheme
               ptn_ad(ji,jj,jk) = ptn_ad(ji,jj,jk) + zhw 
               ptn_ad(ji,jj,jk-1) = ptn_ad(ji,jj,jk-1) + zhw
            END DO
         END DO
      END DO
!$omp end do
!$omp end parallel

      ptn_ad(:,:,1)  = ptn_ad(:,:,1)  + zwzad(:,:,1) * pwn(:,:,1)
      zwzad(:,:,1) = 0._R8

      CALL MYFRTPROF_WALL('TRA_ADV_CEN2V_AD_3D: AD VERT ADVECTION',1)

   END SUBROUTINE tra_adv_cen2_ad_3d

END MODULE traadv_cen2v
