! Call as ZRHO = RHO_UNESCOTL(SB,TB,S,T,0._R8,.FALSE.)
! INTERNAL VARIABLE DECLARATIONS:
!       REAL(R8) ::  A2,A2_TL,B2,B2_TL,C2,ROOTS,RHOW_TL

        RHOW_TL = ( 5._R8*6.536332D-09 * TB**4 -4._R8*1.120083D-06*TB**3 + &
                & 3._R8*1.001685D-04*TB**2 - 2._R8*9.095290D-03*TB + &
                & 6.793952D-02 )
        A2 = (((5.3875D-09 * TB - 8.2467D-07) * TB + 7.6438D-05) &
        &    * TB - 4.0899D-03) * TB + 8.24493D-01
        A2_TL = (4._R8*5.3875D-09 * TB**3 - 3._R8*8.2467D-07 * TB**2 + &
              &  2._R8*7.6438D-05*TB - 4.0899D-03)
        B2 = (-1.6546D-06 * TB + 1.0227D-04) * TB - 5.72466D-03
        B2_TL = (-2._R8*1.6546D-06 * TB + 1.0227D-04)
        C2 = 4.8314D-04

        ROOTS = SQRT (MAX(0._R8,SB))
        ZRHO = 2._R8*C2*SB*S + B2_TL*ROOTS*SB*T + &
        & 1.5_R8*B2*ROOTS*S + A2_TL*SB*T + A2*S + &
        & RHOW_TL*T
