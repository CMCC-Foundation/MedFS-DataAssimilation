        IF( .NOT. LL_SLAAD_INIT ) THEN
           ROOTS = SQRT (MAX(0._R8,SB))
           RHOW_TL = ( 5._R8*6.536332E-09 * TB**4 -4._R8*1.120083E-06*TB**3 + &
                   & 3._R8*1.001685E-04*TB**2 - 2._R8*9.095290E-03*TB + &
                   & 6.793952E-02 )
           A2 = (((5.3875E-09 * TB - 8.2467E-07) * TB + 7.6438E-05) &
           &    * TB - 4.0899E-03) * TB + 8.24493E-01
           A2_TL = (4._R8*5.3875E-09 * TB**3 - 3._R8*8.2467E-07 * TB**2 + &
                 &  2._R8*7.6438E-05*TB - 4.0899E-03)
           B2 = (-1.6546E-06 * TB + 1.0227E-04) * TB - 5.72466E-03
           B2_TL = (-2._R8*1.6546E-06 * TB + 1.0227E-04)
           C2 = 4.8314E-04
           SAUX_AD(K,JVL) = (2._R8*C2*SB + 1.5_R8*B2*ROOTS + A2 )
           TAUX_AD(K,JVL) = (B2_TL*ROOTS*SB + A2_TL*SB + RHOW_TL )
        ENDIF
        S_AD = SAUX_AD(K,JVL)*OBS%GRA(OBS_K)
        T_AD = TAUX_AD(K,JVL)*OBS%GRA(OBS_K)
