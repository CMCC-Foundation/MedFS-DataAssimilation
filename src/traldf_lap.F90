MODULE traldf_lap

   USE TLAD_VARS
   USE MYFRTPROF
   IMPLICIT NONE

  INTERFACE tra_ldf_lap
         MODULE PROCEDURE  tra_ldf_lap_2d, tra_ldf_lap_3d
  END INTERFACE

  INTERFACE tra_ldf_lap_ad
         MODULE PROCEDURE  tra_ldf_lap_ad_2d,tra_ldf_lap_ad_3d
  END INTERFACE

CONTAINS

   SUBROUTINE tra_ldf_lap_3d(jpi,jpj,jpk,ptb,pta) 
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT( in ) ::   jpi,jpj,jpk
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta        ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zabe1, zabe2, zbtr   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ztu, ztv
      !!----------------------------------------------------------------------
      !
      !      e1ur(:,:) = e2u(:,:) / e1u(:,:)
      !      e2vr(:,:) = e1v(:,:) / e2v(:,:)
      !                                                       ! =========== !    
      CALL MYFRTPROF_WALL('tra_ldf_lap_3d: TL DIFFUSION',0)

!$omp parallel  default(shared) private(jk,jj,ji,zabe1,zabe2,zbtr)
!$omp do schedule(dynamic,1)
      DO jk = 1, jpkm1                                            ! slab loop
         !                                           
         ! 1. First derivative (gradient)
         ! -------------------
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               zabe1 = ahtu(ji,jj) * umask(ji,jj,jk) * e1ur(ji,jj) * e3u(jk)
               zabe2 = ahtv(ji,jj) * vmask(ji,jj,jk) * e2vr(ji,jj) * e3v(jk)
               ztu(ji,jj,jk) = zabe1 * ( ptb(ji+1,jj  ,jk) - ptb(ji,jj,jk) )
               ztv(ji,jj,jk) = zabe2 * ( ptb(ji  ,jj+1,jk) - ptb(ji,jj,jk) )
            END DO
         END DO
      
         ! 2. Second derivative (divergence) added to the general tracer trends
         ! ---------------------------------------------------------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zbtr = 1._wp / ( e1t(ji,jj) *e2t(ji,jj) * e3t(jk) )
               ! horizontal diffusive trends added to the general tracer trends
               pta(ji,jj,jk) = pta(ji,jj,jk) + zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk)   &
                  &                                          + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
            END DO
         END DO
         !
      END DO                                             !  End of slab  
!$omp end do
!$omp end parallel
      !
      CALL MYFRTPROF_WALL('tra_ldf_lap_3d: TL DIFFUSION',1)
   END SUBROUTINE tra_ldf_lap_3d

!======================================================================

   SUBROUTINE tra_ldf_lap_2d(jpi,jpj,ptb,pta) 
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT( in ) ::   jpi,jpj
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pta        ! tracer trend 
      !
      INTEGER  ::   ji, jj
      REAL(wp) ::   zabe1, zabe2, zbtr   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   ztu, ztv
      !!----------------------------------------------------------------------
      !
      !      e1ur(:,:) = e2u(:,:) / e1u(:,:)
      !      e2vr(:,:) = e1v(:,:) / e2v(:,:)
      !                                                       ! =========== !    
      !                                           
      ! 1. First derivative (gradient)
      ! -------------------
      DO jj = 1, jpjm1
         DO ji = 1, jpim1   ! vector opt.
            zabe1 = ahtu(ji,jj) * umask(ji,jj,1) * e1ur(ji,jj)
            zabe2 = ahtv(ji,jj) * vmask(ji,jj,1) * e2vr(ji,jj)
            ztu(ji,jj) = zabe1 * ( ptb(ji+1,jj  ) - ptb(ji,jj) )
            ztv(ji,jj) = zabe2 * ( ptb(ji  ,jj+1) - ptb(ji,jj) )
         END DO
      END DO
      
      ! 2. Second derivative (divergence) added to the general tracer trends
      ! ---------------------------------------------------------------------
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            zbtr = 1._wp / ( e1t(ji,jj) *e2t(ji,jj) )
            ! horizontal diffusive trends added to the general tracer trends
            pta(ji,jj) = pta(ji,jj) + zbtr * (  ztu(ji,jj) - ztu(ji-1,jj)   &
               &                                          + ztv(ji,jj) - ztv(ji,jj-1)  )
         END DO
      END DO
      !
   END SUBROUTINE tra_ldf_lap_2d

!======================================================================

   SUBROUTINE tra_ldf_lap_ad_3d(jpi,jpj,jpk,ptb_ad,pta_ad) 
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT( in ) ::   jpi,jpj,jpk
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   ptb_ad     ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta_ad     ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zabe1, zabe2, zbtr, ztaad
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztvad, ztuad
      !!----------------------------------------------------------------------
      !
      !      e1ur(:,:) = e2u(:,:) / e1u(:,:)
      !      e2vr(:,:) = e1v(:,:) / e2v(:,:)
      !                                                       ! =========== !    
      CALL MYFRTPROF_WALL('tra_ldf_lap_AD_3d: AD DIFFUSION',0)
      ztvad(:,:,:) = 0.0_wp
      ztuad(:,:,:) = 0.0_wp
      ztaad        = 0.0_wp
      !
!$omp parallel  default(shared) private(jk,jj,ji,zbtr,ztaad,zabe1,zabe2)
!$omp do schedule(dynamic,1)
      DO jk = 1, jpkm1                                            ! slab loop
         !                                           
         ! 2. Second derivative (divergence) added to the general tracer trends
         ! ---------------------------------------------------------------------
         DO jj = jpjm1,2,-1
            DO ji = jpim1,2,-1
               zbtr = 1._wp / ( e1t(ji,jj) *e2t(ji,jj) * e3t(jk) )
               ztaad= zbtr * pta_ad(ji,jj,jk)
               ! horizontal diffusive trends
               ztuad(ji  ,jj  ,jk) = ztuad(ji  ,jj  ,jk) + ztaad
               ztuad(ji-1,jj  ,jk) = ztuad(ji-1,jj  ,jk) - ztaad
               ztvad(ji  ,jj  ,jk) = ztvad(ji  ,jj  ,jk) + ztaad
               ztvad(ji  ,jj-1,jk) = ztvad(ji  ,jj-1,jk) - ztaad
            END DO
         END DO
         ! 1. First derivative (gradient)
         ! -------------------
         DO jj = jpjm1,1,-1
            DO ji = jpim1,1,-1
               zabe1 = ahtu(ji,jj) * umask(ji,jj,jk) * e1ur(ji,jj) * e3u(jk)
               zabe2 = ahtv(ji,jj) * vmask(ji,jj,jk) * e2vr(ji,jj) * e3v(jk)
               ptb_ad(ji  ,jj  ,jk) = ptb_ad(ji  ,jj  ,jk) - (zabe1 * ztuad(ji,jj,jk) + zabe2 * ztvad(ji,jj,jk))
               ptb_ad(ji+1,jj  ,jk) = ptb_ad(ji+1,jj  ,jk) +  zabe1 * ztuad(ji,jj,jk)
               ptb_ad(ji  ,jj+1,jk) = ptb_ad(ji  ,jj+1,jk) +  zabe2 * ztvad(ji,jj,jk)
               ztuad(ji,jj,jk)   = 0.0_wp
               ztvad(ji,jj,jk)   = 0.0_wp
            END DO
         END DO
        !
      END DO                                             !  End of slab  
!$omp end do
!$omp end parallel
      CALL MYFRTPROF_WALL('tra_ldf_lap_AD_3d: AD DIFFUSION',1)
      !
   END SUBROUTINE tra_ldf_lap_ad_3d

   SUBROUTINE tra_ldf_lap_ad_2d(jpi,jpj,ptb_ad,pta_ad) 
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT( in ) ::   jpi,jpj
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   ptb_ad     ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pta_ad     ! tracer trend 
      !
      INTEGER  ::   ji, jj           ! dummy loop indices
      REAL(wp) ::   zabe1, zabe2, zbtr, ztaad
      REAL(wp), DIMENSION(jpi,jpj) :: ztvad, ztuad
      !!----------------------------------------------------------------------
      !
      !      e1ur(:,:) = e2u(:,:) / e1u(:,:)
      !      e2vr(:,:) = e1v(:,:) / e2v(:,:)
      !                                                       ! =========== !    
      ztvad(:,:) = 0.0_wp
      ztuad(:,:) = 0.0_wp
      ztaad        = 0.0_wp
      !
      !                                           
      ! 2. Second derivative (divergence) added to the general tracer trends
      ! ---------------------------------------------------------------------
      DO jj = jpjm1,2,-1
         DO ji = jpim1,2,-1
            zbtr = 1._wp / ( e1t(ji,jj) *e2t(ji,jj) )
            ztaad= zbtr * pta_ad(ji,jj)
            ! horizontal diffusive trends
            ztuad(ji  ,jj  ) = ztuad(ji  ,jj  ) + ztaad
            ztuad(ji-1,jj  ) = ztuad(ji-1,jj  ) - ztaad
            ztvad(ji  ,jj  ) = ztvad(ji  ,jj  ) + ztaad
            ztvad(ji  ,jj-1) = ztvad(ji  ,jj-1) - ztaad
         END DO
      END DO
      !
      ! 1. First derivative (gradient)
      ! -------------------
      DO jj = jpjm1,1,-1
         DO ji = jpim1,1,-1
            zabe1 = ahtu(ji,jj) * umask(ji,jj,1) * e1ur(ji,jj) 
            zabe2 = ahtv(ji,jj) * vmask(ji,jj,1) * e2vr(ji,jj)
            ptb_ad(ji  ,jj  ) = ptb_ad(ji  ,jj  ) - (zabe1 * ztuad(ji,jj) + zabe2 * ztvad(ji,jj))
            ptb_ad(ji+1,jj  ) = ptb_ad(ji+1,jj  ) +  zabe1 * ztuad(ji,jj)
            ptb_ad(ji  ,jj+1) = ptb_ad(ji  ,jj+1) +  zabe2 * ztvad(ji,jj)
            ztuad(ji,jj)   = 0.0_wp
            ztvad(ji,jj)   = 0.0_wp
         END DO
      END DO
      !
   END SUBROUTINE tra_ldf_lap_ad_2d

END MODULE traldf_lap
