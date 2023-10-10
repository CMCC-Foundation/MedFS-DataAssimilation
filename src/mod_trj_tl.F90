      SUBROUTINE MOD_TRJ_TL( JPI,JPJ,UMOD,VMOD,E1U,E2V,              &
                             JPT,JPN,PIMOD,PJMOD,PTIME,UTAN,VTAN,XTAN,YTAN )
!!----------------------------------------------------------------------
!! ARGUMENTS
!! =========
      IMPLICIT NONE
      INTEGER*8 :: JPT,JPN
      INTEGER :: JPI,JPJ
      REAL*8, DIMENSION(JPT+1,JPN) :: PIMOD,PJMOD
      REAL*8, DIMENSION(JPN) :: PTIME
      REAL*8, DIMENSION(JPI,JPJ) :: UMOD,VMOD,E1U,E2V
      REAL*8, DIMENSION(JPN) :: XTAN,YTAN
!!----------------------------------------------------------------------
!! LOCAL DECLARATIONS
!! ==================
      INTEGER :: JN,JT
      REAL*8, DIMENSION(JPI,JPJ) :: UTAN,VTAN
      REAL*8 :: ZPIN,ZPJN,ZPITAN,ZPJTAN,ZRDT,ZUTAN,ZVTAN
!!----------------------------------------------------------------------

!*... TRAJECTORIES COMPUTATION
      DO JN = 1, JPN
        ZPITAN = 0.E0
        ZPJTAN = 0.E0
        ZRDT = PTIME(JN)*3600./FLOAT(JPT)
        DO JT = 1, JPT
          ZPIN = PIMOD(JT,JN)
          ZPJN = PJMOD(JT,JN)
          CALL FLOITPTAN( JPI,JPJ,UMOD,VMOD,UTAN,VTAN,E1U,E2V,    &
                          ZPIN,ZPJN,ZPITAN,ZPJTAN,ZUTAN,ZVTAN )
          ZUTAN  = ZPITAN + ZRDT * ZUTAN
          ZVTAN  = ZPJTAN + ZRDT * ZVTAN
          ZPITAN = ZUTAN
          ZPJTAN = ZVTAN
        END DO
        XTAN(JN) = ZPITAN
        YTAN(JN) = ZPJTAN
      END DO

      RETURN
      END SUBROUTINE MOD_TRJ_TL

!!======================================================================

      SUBROUTINE FLOITPTAN( JPI,JPJ,UMOD,VMOD,UTAN,VTAN,E1U,E2V,     &
                            PIFL,PJFL,PIFLTL,PJFLTL,PUFLTL,PVFLTL )
!!----------------------------------------------------------------------
!! ARGUMENTS
!! =========
      IMPLICIT NONE
      INTEGER :: JPI,JPJ
      REAL*8 :: PIFL,PJFL,PIFLTL,PJFLTL,PUFLTL,PVFLTL
      REAL*8, DIMENSION(JPI,JPJ) :: UMOD,VMOD,UTAN,VTAN,E1U,E2V
!!----------------------------------------------------------------------
!! LOCAL DECLARATIONS
!! ==================
      INTEGER :: IIL,IJL,JIND1,JIND2
      INTEGER, DIMENSION(2) :: IID,IJD
      REAL*8, DIMENSION(2) :: ZLAGX,ZLAGY,ZLAGXTL,ZLAGYTL
      REAL*8, DIMENSION(2,2) :: ZUV,ZUVTL
!!----------------------------------------------------------------------

!! 1. INTERPOLATION OF THE ZONAL VELOCITY
!! ======================================

!*... NEIGHBOORING POINTS (BACKGROUND)
      IIL = INT(PIFL-.5)
      IJL = INT(PJFL   )
      DO JIND1 = 1, 2
        IID(JIND1) = IIL + JIND1 - 1
        IJD(JIND1) = IJL + JIND1 - 1
      END DO

!*... LAGRANGE COEFFICIENTS (BACKGROUND)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        IF( JIND1.NE.JIND2 ) THEN
          ZLAGX(JIND1) = ( PIFL - (FLOAT(IID(JIND2))+.5) ) / FLOAT( IID(JIND1)-IID(JIND2) )
          ZLAGY(JIND1) = ( PJFL -  FLOAT(IJD(JIND2))     ) / FLOAT( IJD(JIND1)-IJD(JIND2) )
        ENDIF
      END DO
      END DO

!*... VALUE OF THE ZONAL VELOCITY (BACKGROUND)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        ZUV(JIND1,JIND2) = UMOD(IID(JIND1),IJD(JIND2)) / E1U(IID(JIND1),IJD(JIND2))
      END DO
      END DO

!*... LAGRANGE COEFFICIENTS (TANGENT)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        IF( JIND1.NE.JIND2 ) THEN
          ZLAGXTL(JIND1) = PIFLTL / FLOAT( IID(JIND1)-IID(JIND2) )
          ZLAGYTL(JIND1) = PJFLTL / FLOAT( IJD(JIND1)-IJD(JIND2) )
        ENDIF
      END DO
      END DO

!*... VALUE OF THE ZONAL VELOCITY (TANGENT)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        ZUVTL(JIND1,JIND2) = UTAN(IID(JIND1),IJD(JIND2)) / E1U(IID(JIND1),IJD(JIND2))
      END DO
      END DO

!*... INTERPOLATION OF THE ZONAL VELOCITY
      PUFLTL = 0.E0
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        PUFLTL = PUFLTL                                                  &
               + ZUV(JIND1,JIND2)   * ZLAGXTL(JIND1) * ZLAGY(JIND2)      &
               + ZUV(JIND1,JIND2)   * ZLAGX(JIND1)   * ZLAGYTL(JIND2)    &
               + ZUVTL(JIND1,JIND2) * ZLAGX(JIND1)   * ZLAGY(JIND2)
      END DO
      END DO

!! 2. INTERPOLATION OF THE MERIDIAN VELOCITY
!! =========================================

!*... NEIGHBOORING POINTS (BACKGROUND)
      IIL = INT(PIFL   )
      IJL = INT(PJFL-.5)
      DO JIND1 = 1, 2
        IID(JIND1) = IIL + JIND1 - 1
        IJD(JIND1) = IJL + JIND1 - 1
      END DO

!*... LAGRANGE COEFFICIENTS (BACKGROUND)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        IF( JIND1.NE.JIND2 ) THEN
          ZLAGX(JIND1) = ( PIFL -  FLOAT(IID(JIND2))     ) / FLOAT( IID(JIND1)-IID(JIND2) )
          ZLAGY(JIND1) = ( PJFL - (FLOAT(IJD(JIND2))+.5) ) / FLOAT( IJD(JIND1)-IJD(JIND2) )
        ENDIF
      END DO
      END DO

!*... VALUE OF THE MERIDIAN VELOCITY (BACKGROUND)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        ZUV(JIND1,JIND2) = VMOD(IID(JIND1),IJD(JIND2)) / E2V(IID(JIND1),IJD(JIND2))
      END DO
      END DO

!*... LAGRANGE COEFFICIENTS (TANGENT)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        IF( JIND1.NE.JIND2 ) THEN
          ZLAGXTL(JIND1) = PIFLTL / FLOAT( IID(JIND1)-IID(JIND2) )
          ZLAGYTL(JIND1) = PJFLTL / FLOAT( IJD(JIND1)-IJD(JIND2) )
        ENDIF
      END DO
      END DO

!*... VALUE OF THE MERIDIAN VELOCITY (TANGENT)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        ZUVTL(JIND1,JIND2) = VTAN(IID(JIND1),IJD(JIND2)) / E2V(IID(JIND1),IJD(JIND2))
      END DO
      END DO

!*... INTERPOLATION OF THE MERIDIAN VELOCITY
      PVFLTL = 0.E0
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        PVFLTL = PVFLTL                                                 &
               + ZUV(JIND1,JIND2)   * ZLAGXTL(JIND1) * ZLAGY(JIND2)     &
               + ZUV(JIND1,JIND2)   * ZLAGX(JIND1)   * ZLAGYTL(JIND2)   &
               + ZUVTL(JIND1,JIND2) * ZLAGX(JIND1)   * ZLAGY(JIND2)
      END DO
      END DO

      RETURN
      END SUBROUTINE FLOITPTAN
