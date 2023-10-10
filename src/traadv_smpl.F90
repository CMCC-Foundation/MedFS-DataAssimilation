MODULE traadv_smpl

   USE TLAD_VARS
   USE MYFRTPROF
   IMPLICIT NONE

   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zwx, zwy, zwz
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwxad, zwyad, zwzad
   LOGICAL, PRIVATE, SAVE :: LLINIT_TL = .FALSE.
   LOGICAL, PRIVATE, SAVE :: LLINIT_AD = .FALSE.

  INTERFACE tra_adv_cen2
         MODULE PROCEDURE  tra_adv_cen2_2d, tra_adv_cen2_3d
  END INTERFACE

  INTERFACE tra_adv_cen2_ad
         MODULE PROCEDURE  tra_adv_cen2_ad_2d, tra_adv_cen2_ad_3d
  END INTERFACE

CONTAINS

   SUBROUTINE tra_adv_smpl(jpi,jpj,jpk,pun,pvn,pwn,ptn,pta)

      IMPLICIT NONE

      !
      INTEGER, INTENT( in ) ::   jpi,jpj,jpk
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: pun,pvn,pwn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: ptn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(INOUT) :: pta

      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('tra_adv_smpl: TL ADVECTION',0)
      IF( .NOT. LLINIT_TL ) THEN
        ALLOCATE( zwx(jpi,jpj), &
                zwy(jpi,jpj), &
                zwz(jpi,jpj)  )
        LLINIT_TL = .TRUE.
      ENDIF
      zwx = 0._wp
      zwy = 0._wp
      zwz = 0._wp

      DO jk=1,km
         DO jj=1,jpj
            DO ji=1,jpi
               zwx(ji,jj) = pun(ji,jj,jk)*ptn(i+1,j,k)*e2t(i,j)
            ENDDO
         ENDDO

         DO jj=1,jpj
            DO ji=1,jpi
               zwy(ji,jj) = pvn(ji,jj,jk)*ptn(i,j+1,k)*e1t(i,j)
            ENDDO
         ENDDO

         DO jj=1,jpj
            DO ji=1,jpi
               pta(ji,jj,jk) = pta(ji,jj,jk) - 
               ( zwx(ji,jj) - zwx(ji-1,jj) + zwy(ji,jj) - zwy(ji,jj-1) )       &
                                          / ( e1t(ji,jj) * e2t(ji,jj) )
            ENDDO
         ENDDO
      ENDDO
      !                                          ! set time step
      CALL MYFRTPROF_WALL('tra_adv_smpl: TL ADVECTION',1)

   END SUBROUTINE tra_adv_smpl

   SUBROUTINE tra_adv_smplv(jpi,jpj,jpk,pun,pvn,pwn,ptn,pta)

      IMPLICIT NONE

      !
      INTEGER, INTENT( in ) ::   jpi,jpj,jpk
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: pun,pvn,pwn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: ptn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(INOUT) :: pta

      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('tra_adv_smplv: TL ADVECTION',0)
      IF( .NOT. LLINIT_TL ) THEN
        ALLOCATE( zwz(jpi,jpj,jpk) )
        LLINIT_TL = .TRUE.
      ENDIF
      zwz = 0._wp

      DO jk=1,km
         DO jj=1,jpj
            DO ji=1,jpi
               zwx(ji,jj) = pun(ji,jj,jk)*ptn(i+1,j,k)*e2t(i,j)
            ENDDO
         ENDDO

         DO jj=1,jpj
            DO ji=1,jpi
               zwy(ji,jj) = pvn(ji,jj,jk)*ptn(i,j+1,k)*e1t(i,j)
            ENDDO
         ENDDO

         DO jj=1,jpj
            DO ji=1,jpi
               pta(ji,jj,jk) = pta(ji,jj,jk) -
               ( zwx(ji,jj) - zwx(ji-1,jj) + zwy(ji,jj) - zwy(ji,jj-1) )       &
                                          / ( e1t(ji,jj) * e2t(ji,jj) )
            ENDDO
         ENDDO
      ENDDO
      !                                          ! set time step
      CALL MYFRTPROF_WALL('tra_adv_smplv: TL ADVECTION',1)

   END SUBROUTINE tra_adv_smplv


   SUBROUTINE tra_adv_smpl_ad(jpi,jpj,jpk,pun,pvn,pwn,ptn_ad,pta_ad)

      IMPLICIT NONE

      !
      INTEGER, INTENT( in ) ::   jpi,jpj,jpk
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: pun,pvn,pwn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(INOUT) :: ptn_ad
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(INOUT) :: pta_ad

      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   ierr             ! local integer
      REAL(wp) ::   zbtr, ztra                            ! local scalars
      REAL(wp) ::   zfui, zfvj, zhw
      !!----------------------------------------------------------------------
      !
      CALL MYFRTPROF_WALL('tra_adv_cen2_AD_3d: AD ADVECTION',0)

      !
      IF( .NOT. LLINIT_AD ) THEN
        ALLOCATE( zwxad(jpi,jpj,jpk), &
                zwyad(jpi,jpj,jpk), &
                zwzad(jpi,jpj,jpk)  )
        LLINIT_AD = .TRUE.
      ENDIF
      zwxad = 0._R8
      zwyad = 0._R8
      zwzad = 0._R8
      !                                          ! set time step
      ! II. Divergence of advective fluxes
      ! ----------------------------------
!__NOT-USED__!  parallel  default(shared) private(jk,jj,ji,zbtr)
!__NOT-USED__!  do schedule(dynamic,1)
      DO jk = jpkm1, 1, -1
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) *  e3t(jk) )
               zwxad(ji,jj,jk) = zwxad(ji,jj,jk) - zbtr * pta_ad(ji,jj,jk)
               zwyad(ji,jj,jk) = zwyad(ji,jj,jk) - zbtr * pta_ad(ji,jj,jk)
               zwzad(ji,jj,jk) = zwzad(ji,jj,jk) - zbtr * pta_ad(ji,jj,jk)
               zwxad(ji-1,jj,jk) = zwxad(ji-1,jj,jk) + zbtr * pta_ad(ji,jj,jk)
               zwyad(ji,jj-1,jk) = zwyad(ji,jj-1,jk) + zbtr * pta_ad(ji,jj,jk)
               zwzad(ji,jj,jk+1) = zwzad(ji,jj,jk+1) + zbtr * pta_ad(ji,jj,jk)
            END DO
         END DO
      END DO
!__NOT-USED__!  end do
!__NOT-USED__!  end parallel

!__NOT-USED__!  parallel  default(shared) private(jj,jk,ji,zhw,zbtr)
!__NOT-USED__!  do schedule(dynamic,1)
      !!!!!... ptn_ad = 0._R8
      DO jk = jpk, 2, -1      ! Second order centered tracer flux at w-point
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               zhw   = 0.5 * pwn(ji,jj,jk)*zwzad(ji,jj,jk)
               zwzad(ji,jj,jk) = 0._R8
               ! centered scheme
               ptn_ad(ji,jj,jk) = ptn_ad(ji,jj,jk) + zhw 
               ptn_ad(ji,jj,jk-1) = ptn_ad(ji,jj,jk-1) + zhw
            END DO
         END DO
      END DO
!__NOT-USED__!  end do
!__NOT-USED__!  end parallel

      ptn_ad(:,:,1)  = ptn_ad(:,:,1)  + zwzad(:,:,1) * pwn(:,:,1)

!__NOT-USED__!  parallel  default(shared) private(jk,jj,ji,zfui,zfvj)
!__NOT-USED__!  do schedule(dynamic,1)
      DO jk = jpkm1, 1, -1
         !                        ! Second order centered tracer flux at u- and v-points
         DO jj = jpjm1, 1, -1
            DO ji = jpim1, 1, -1   ! vector opt.
               ! volume fluxes * 1/2
               zfvj = 0.5 * pvn(ji,jj,jk)*zwyad(ji,jj,jk)
               zwyad(ji,jj,jk) = 0._wp
               zfui = 0.5 * pun(ji,jj,jk)*zwxad(ji,jj,jk)
               zwxad(ji,jj,jk) = 0._wp
               ! centered scheme
               ptn_ad(ji  ,jj  ,jk) = ptn_ad(ji  ,jj  ,jk) + zfvj
               ptn_ad(ji  ,jj+1,jk) = ptn_ad(ji  ,jj+1,jk) + zfvj
               ptn_ad(ji  ,jj  ,jk) = ptn_ad(ji  ,jj  ,jk) + zfui
               ptn_ad(ji+1,jj  ,jk) = ptn_ad(ji+1,jj  ,jk) + zfui
            END DO
         END DO
      END DO
!__NOT-USED__!  end do
!__NOT-USED__!  end parallel

      CALL MYFRTPROF_WALL('tra_adv_smpl_AD: AD ADVECTION',1)

   END SUBROUTINE tra_adv_smpl_ad

END MODULE traadv_smpl
