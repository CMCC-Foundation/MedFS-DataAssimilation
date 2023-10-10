      SUBROUTINE MOD_TRJ_AD( JPI,JPJ,UMOD,VMOD,E1U,E2V,        &
                             JPT,JPN,PIMOD,PJMOD,PTIME,UADJ,VADJ,XADJ,YADJ )
!!----------------------------------------------------------------------
!! ARGUMENTS
!! =========
      IMPLICIT NONE
      INTEGER*8 :: JPT,JPN
      INTEGER :: JPI,JPJ
      REAL*8, DIMENSION(JPT+1,JPN) :: PIMOD,PJMOD
      REAL*8, DIMENSION(JPN) :: PTIME
      REAL*8, DIMENSION(JPI,JPJ) :: UMOD,VMOD,E1U,E2V
      REAL*8, DIMENSION(JPN) :: XADJ,YADJ
!!----------------------------------------------------------------------
!! LOCAL DECLARATIONS
!! ==================
      INTEGER :: JN,JT
      REAL*8, DIMENSION(JPI,JPJ) :: UADJ,VADJ
      REAL*8 :: ZPIN,ZPJN,ZUADJ,ZVADJ,ZRDT,ZPIADJ,ZPJADJ
!!----------------------------------------------------------------------

!*... TRAJECTORIES COMPUTATION
      UADJ(:,:) = 0.E0
      VADJ(:,:) = 0.E0
      DO JN = 1, JPN
        ZPIADJ = XADJ(JN)
        ZPJADJ = YADJ(JN)
        ZRDT = PTIME(JN)*3600./FLOAT(JPT)
        DO JT = JPT, 1, -1
          ZPIN = PIMOD(JT,JN)
          ZPJN = PJMOD(JT,JN)
          ZUADJ  = ZPIADJ
          ZVADJ  = ZPJADJ
          ZUADJ  = ZUADJ * ZRDT
          ZVADJ  = ZVADJ * ZRDT
          CALL FLOITPADJ( JPI,JPJ,UMOD,VMOD,UADJ,VADJ,E1U,E2V, &
                          ZPIN,ZPJN,ZPIADJ,ZPJADJ,ZUADJ,ZVADJ )
        END DO
      END DO

      RETURN
      END SUBROUTINE MOD_TRJ_AD

!!======================================================================

      SUBROUTINE FLOITPADJ( JPI,JPJ,UMOD,VMOD,UADJ,VADJ,E1U,E2V, &
                             PIFL,PJFL,PIFLAD,PJFLAD,PUFLAD,PVFLAD )
!!----------------------------------------------------------------------
!! ARGUMENTS
!! =========
      IMPLICIT NONE
      INTEGER :: JPI,JPJ
      REAL*8 :: PIFL,PJFL,PIFLAD,PJFLAD,PUFLAD,PVFLAD
      REAL*8, DIMENSION(JPI,JPJ) :: UMOD,VMOD,UADJ,VADJ,E1U,E2V
!!----------------------------------------------------------------------
!! LOCAL DECLARATIONS
!! ==================
      INTEGER :: IIL,IJL,JIND1,JIND2
      INTEGER, DIMENSION(2) :: IID,IJD
      REAL*8, DIMENSION(2) :: ZLAGX,ZLAGY,ZLAGXAD,ZLAGYAD
      REAL*8, DIMENSION(2,2) :: ZUV,ZUVAD
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

!*... INTERPOLATION OF THE ZONAL VELOCITY
      ZLAGXAD(:) = 0.E0
      ZLAGYAD(:) = 0.E0
      ZUVAD(:,:) = 0.E0
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        ZLAGXAD(JIND1)     = ZLAGXAD(JIND1) + PUFLAD * ZUV(JIND1,JIND2) * ZLAGY(JIND2)
        ZLAGYAD(JIND2)     = ZLAGYAD(JIND2) + PUFLAD * ZUV(JIND1,JIND2) * ZLAGX(JIND1)
        ZUVAD(JIND1,JIND2) = ZUVAD(JIND1,JIND2) + PUFLAD * ZLAGX(JIND1) * ZLAGY(JIND2)
      END DO
      END DO
      PUFLAD = 0.E0

!*... VALUE OF THE ZONAL VELOCITY (ADJOINT)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        UADJ(IID(JIND1),IJD(JIND2)) = UADJ(IID(JIND1),IJD(JIND2)) &
                                    + ZUVAD(JIND1,JIND2) / E1U(IID(JIND1),IJD(JIND2))
      END DO
      END DO

!*... LAGRANGE COEFFICIENTS (ADJOINT)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        IF( JIND1.NE.JIND2 ) THEN
          PIFLAD = PIFLAD + ZLAGXAD(JIND1) / FLOAT( IID(JIND1)-IID(JIND2) )
          PJFLAD = PJFLAD + ZLAGYAD(JIND1) / FLOAT( IID(JIND1)-IID(JIND2) )
        ENDIF
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

!*... INTERPOLATION OF THE MERIDIAN VELOCITY
      ZLAGXAD(:) = 0.E0
      ZLAGYAD(:) = 0.E0
      ZUVAD(:,:) = 0.E0
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        ZLAGXAD(JIND1)     = ZLAGXAD(JIND1) + PVFLAD * ZUV(JIND1,JIND2) * ZLAGY(JIND2)
        ZLAGYAD(JIND2)     = ZLAGYAD(JIND2) + PVFLAD * ZUV(JIND1,JIND2) * ZLAGX(JIND1)
        ZUVAD(JIND1,JIND2) = ZUVAD(JIND1,JIND2) + PVFLAD * ZLAGX(JIND1) * ZLAGY(JIND2)
      END DO
      END DO
      PVFLAD = 0.E0

!*... VALUE OF THE MERIDIAN VELOCITY (ADJOINT)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        VADJ(IID(JIND1),IJD(JIND2)) = VADJ(IID(JIND1),IJD(JIND2)) &
                                    + ZUVAD(JIND1,JIND2) / E2V(IID(JIND1),IJD(JIND2))
      END DO
      END DO

!*... LAGRANGE COEFFICIENTS (ADJOINT)
      DO JIND1 = 1, 2
      DO JIND2 = 1, 2
        IF( JIND1.NE.JIND2 ) THEN
          PIFLAD = PIFLAD + ZLAGXAD(JIND1) / FLOAT( IID(JIND1)-IID(JIND2) )
          PJFLAD = PJFLAD + ZLAGYAD(JIND1) / FLOAT( IID(JIND1)-IID(JIND2) )
        ENDIF
      END DO
      END DO

      RETURN
      END SUBROUTINE FLOITPADJ
