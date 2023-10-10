   SUBROUTINE COFDIS_2D(JPI, JPJ, &
   & GLAMT, GPHIT, &
   & TMASK, PDCT )

   USE SET_KND

   ! PROVIDES 2D DISTANCE IN KM
   ! 
   ! A.S. ADAPTED FROM NEMO 3.3
   ! SEMPLIFIED VERSION FOR NPERIO==3
   ! AND CONSIDERING ONLY T-POINTS

   IMPLICIT NONE

      INTEGER(I4) , INTENT( IN)  ::   JPI, JPJ
      REAL(R8), DIMENSION(JPI,JPJ), INTENT( IN)  ::   GLAMT, GPHIT, &
      & TMASK
      REAL(R8), DIMENSION(JPI,JPJ), INTENT( OUT ) ::   PDCT
      !!
      INTEGER(I4) ::   JI, JJ, JL      ! DUMMY LOOP INDICES
      INTEGER(I4) ::   IJU, IJT            ! TEMPORARY INTEGERS
      INTEGER(I4) ::   ICOAST, ITIME
      INTEGER(I4) ::   ICOT, JPIM1, JPJM1
      LOGICAL, DIMENSION(JPI,JPJ) ::   LLCOTF   ! Coastal points
      REAL(R8) ::   ZDATE0
      REAL(R8), PARAMETER :: RA = 6371.229_R8
      REAL(R8), PARAMETER :: RAD = 3.141592653589793_R8 / 180._R8

      ! CARTESIAN COORDINATES FOR T-POINTS
      REAL(R8), DIMENSION(JPI,JPJ)   ::   ZXT, ZYT, ZZT, ZMASK   

      ! TEMPORARY WORKSPACE
      REAL(R8), DIMENSION(3*JPI*JPJ) ::   ZXC, ZYC, ZZC, ZDIS    
      !!----------------------------------------------------------------------

      PDCT(:,:) = 0._R8
      ZXT(:,:) = COS( RAD * GPHIT(:,:) ) * COS( RAD * GLAMT(:,:) )
      ZYT(:,:) = COS( RAD * GPHIT(:,:) ) * SIN( RAD * GLAMT(:,:) )
      ZZT(:,:) = SIN( RAD * GPHIT(:,:) )

      JPIM1 = JPI - 1
      JPJM1 = JPJ - 1

         ! DEFINE THE COASTLINE POINTS (U, V AND F)
         DO JJ = 2, JPJM1
            DO JI = 2, JPIM1
               ZMASK(JI,JJ) =  ( TMASK(JI,JJ+1) + TMASK(JI+1,JJ+1) &
                   &           + TMASK(JI,JJ  ) + TMASK(JI+1,JJ  ) )
               LLCOTF(JI,JJ) = ( ZMASK(JI,JJ)>0._R8 ).AND.(ZMASK(JI,JJ)<4._R8 )
            END DO
         END DO

         ! LATERAL BOUNDARIES CONDITIONS
         LLCOTF(:, 1 ) = TMASK(:,  2  ) == 1
         LLCOTF(:,JPJ) = TMASK(:,JPJM1) == 1

         LLCOTF( 1 ,:) = LLCOTF(JPIM1,:)
         LLCOTF(JPI,:) = LLCOTF(  2  ,:)

            DO JI = 1, JPIM1
               IJU = JPI - JI + 1
               LLCOTF(JI,JPJM1) = LLCOTF(IJU,JPJ-2)
               LLCOTF(JI,JPJ  ) = LLCOTF(IJU,JPJ-3)
            END DO

         ! COMPUTE CARTESIAN COORDINATES OF COASTLINE POINTS
         ! AND THE NUMBER OF COASTLINE POINTS
         ICOAST = 0
         DO JJ = 1, JPJ
            DO JI = 1, JPI
               IF( LLCOTF(JI,JJ) ) THEN
                  ICOAST = ICOAST + 1
                  ZXC(ICOAST) = COS( RAD*GPHIT(JI,JJ) ) * COS( RAD*GLAMT(JI,JJ) )
                  ZYC(ICOAST) = COS( RAD*GPHIT(JI,JJ) ) * SIN( RAD*GLAMT(JI,JJ) )
                  ZZC(ICOAST) = SIN( RAD*GPHIT(JI,JJ) )
               ENDIF
            END DO
         END DO

         ! DISTANCE FOR THE T-POINTS
         DO JJ = 1, JPJ
            DO JI = 1, JPI
               IF( TMASK(JI,JJ) == 0._R8 ) THEN
                  PDCT(JI,JJ) = 0._R8
               ELSE
                  DO JL = 1, ICOAST
                     ZDIS(JL) = ( ZXT(JI,JJ) - ZXC(JL) )**2   &
                        &     + ( ZYT(JI,JJ) - ZYC(JL) )**2   &
                        &     + ( ZZT(JI,JJ) - ZZC(JL) )**2
                  END DO
                  PDCT(JI,JJ) = RA * SQRT( MINVAL( ZDIS(1:ICOAST) ) )
               ENDIF
            END DO
         END DO

   END SUBROUTINE COFDIS_2D
