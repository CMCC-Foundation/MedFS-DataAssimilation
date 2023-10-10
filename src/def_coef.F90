SUBROUTINE DEF_COEF(SIGMA,B1,B2,B3,B)

  USE SET_KND

       IMPLICIT NONE

       REAL(R8)  SIGMA,B0,B1,B2,B3,B,Q

       IF( SIGMA .GE. 2.5_R8 ) THEN
          Q = 0.98711_R8 * SIGMA - 0.96330_R8
       ELSEIF( SIGMA .LT. 2.5_R8 ) THEN
          Q = 3.97156_R8 - 4.14554_R8 *SQRT( 1._R8 - 0.26891_R8* SIGMA )
       ENDIF

       B0 = 1.57825_R8 + (2.44413_R8 * Q) + (1.4281_R8 * Q**2) + (0.422205_R8 * Q**3)
       B1 =( (2.44413_R8 * Q)  + (2.85619_R8 * Q**2) + (1.26661_R8 * Q**3))/B0
       B2 =( -1.4281_R8 * Q**2 -1.26661_R8 * Q**3 )/B0
       B3 =( 0.422205_R8 * Q**3)/B0
       B  = 1._R8 -  ( B1 + B2 + B3 )

END SUBROUTINE DEF_COEF
