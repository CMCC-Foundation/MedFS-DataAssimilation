MODULE DRWN
  !!
  !!
  USE SET_KND
  IMPLICIT NONE
  !!
  !!
  PRIVATE
  !!
  PUBLIC :: DROWN, FILL_EXTRA_BANDS, EXTRA_2_EAST, EXTRA_2_WEST
  !!
  !!
  LOGICAL, PARAMETER :: LSMOOTH = .TRUE.
  !!
CONTAINS
  !!
  !!
  SUBROUTINE DROWN(K_EW, NI, NJ, X, MASK)
    !!
    !!#############################################################################
    !!
    !!  PURPOSE : FILL LAND VALUES ONLY WITH SURROUNDING SEA VALUES
    !!  --------
    !!  NI  : X DIMENSION OF ARRAY X AND MASK            (INTEGER)     
    !!  NJ  : Y DIMENSION OF ARRAY X AND MASK            (INTEGER)
    !!  X   : CURRENTLY TREATED ARRAY                    (2D ARRAY)
    !!  MASK : LAND-SEA MASK                              (2D ARRAY)
    !!
    !!#############################################################################
    !!
    !!
    !! ARGUMENTS :
    !! -----------
    INTEGER(I4),                        INTENT(IN)    :: K_EW, NI, NJ
    !!  
    REAL(R8), DIMENSION(NI,NJ), INTENT(INOUT) :: X
    !!
    INTEGER(I4), DIMENSION(NI,NJ),      INTENT(IN)    :: MASK
    !!
    !!
    !! LOCAL :
    !! --------
    INTEGER(I4), ALLOCATABLE, DIMENSION(:,:) :: MASKV, MASK2 
    !!
    REAL(R8), ALLOCATABLE, DIMENSION(:,:) :: DATA_OLD  
    !!
    INTEGER(I4) :: &
         &      JINC,          &
         &      JI, JJ, JT,    &
         &      NCPT, NSEA,    &
         &      ND_NORTH,      &
         &      ND_EAST,       &
         &      ND_SOUTH,      &
         &      ND_WEST, KK
    !!  
    CHARACTER(LEN=8) :: QUER
    !!
    INTEGER(I4), PARAMETER ::   &
         &      NBAND = 2,  &
         &      NPROB = 1     !: ITERATIVE DISTANCE TO PROBE FOR SEA PRESENC
    !!
    !!
    IF ( SUM(MASK) == 0 ) THEN
       PRINT *, 'THE MASK DOES NOT HAVE SEA POINTS! SKIPPING DROWN!'
       RETURN
    ELSE
    !!
    !! BACKING UP ORIGINAL MASK INTO MASK2(:,:)
    ALLOCATE ( MASKV(NI,NJ), MASK2(NI,NJ), DATA_OLD(NI,NJ) )
    MASKV = MASK
    MASK2 = MASK
    !!
    !!
    JINC = 0
    DO WHILE ( ANY(MASKV(2:NI-2,2:NJ-2) == 0) ) ! AS LONG AS WE STILL HAVE CONTINENTAL POINTS
       JINC = JINC + 1
       !!
       MASKV    = MASK2
       DATA_OLD = X
       !!
       DO JJ=NBAND, NJ-NBAND+1
          JI = 1
          DO WHILE ( JI <= NI )
             IF( MASKV(JI,JJ) == 0 ) THEN
                CALL TEST4SEA(JI,JJ,NI,NJ,DATA_OLD,MASKV,NPROB,1,NSEA)
                ND_NORTH = NSEA
                !!              
                CALL TEST4SEA(JI,JJ,NI,NJ,DATA_OLD,MASKV,NPROB,2,NSEA)
                ND_EAST = NSEA
                !!              
                CALL TEST4SEA(JI,JJ,NI,NJ,DATA_OLD,MASKV,NPROB,3,NSEA)
                ND_SOUTH = NSEA
                !!              
                CALL TEST4SEA(JI,JJ,NI,NJ,DATA_OLD,MASKV,NPROB,4,NSEA)
                ND_WEST = NSEA
                !!              
                IF((ND_NORTH==-1).AND.(ND_EAST==-1).AND.  &
                     &  (ND_SOUTH==-1).AND.(ND_WEST ==-1)) THEN
                   JI = JI+1       ! WE ARE ON EARTH FAR FROM SEA (AT LEAST NPROB POINTS DISTANT)... 
                   ! PUT INITIAL LAND VALUES
                ELSE
                   CALL DECISION(JI,JJ,NI,NJ,X,MASKV,MASK2) 
                 ! WE ARE ON EARTH BUT AT LEAST A SEA POINT IS CLOSER THAN NPROB
                   JI = JI+1
                ENDIF
             ELSE
                JI = JI+1              ! WE ARE ON SEA, JUST MOVING ON...
             ENDIF
          ENDDO
          !!        
       ENDDO
       !!
    ENDDO
    !!
    !! FIRST AND LAST LINE SOMETIMES SHOW WEIRD FEATURES...
    DO JI = 1, NI
       IF ( MASK(JI,1)  == 0 ) X(JI,1)  = X(JI,2)
       IF ( MASK(JI,NJ) == 0 ) X(JI,NJ) = X(JI,NJ-1)
    END DO
    !!
    !! SMOOTHING THE WHAT'S BEEN DONE ON LAND:
    !!
    IF ( LSMOOTH ) THEN
       DO KK = 1, 50
          !!
          DO JJ = 2, NJ-1
             !!
             DO JI = 2, NI-1
                IF ( MASK(JI,JJ) == 0 )   X(JI,JJ) = &
                     &    0.35*X(JI,JJ) + 0.65*0.125*( &
                     &    X(JI+1,JJ)   + X(JI,JJ+1)   + X(JI-1,JJ)   + X(JI,JJ-1) + &
                     &    X(JI+1,JJ+1) + X(JI-1,JJ+1) + X(JI-1,JJ-1) + X(JI+1,JJ-1)   )
             END DO
             !!
             IF (K_EW /= -1) THEN   ! WE CAN USE EAST-WEST PERIODICITY
                IF ( MASK(1,JJ) == 0 )   X(1,JJ) = &
                     &    0.35*X(1,JJ) + 0.65*0.125*( &
                     &    X(2,JJ)   + X(1,JJ+1)   + X(NI-K_EW,JJ)   + X(1,JJ-1) + &
                     &    X(2,JJ+1) + X(NI-K_EW,JJ+1) + X(NI-K_EW,JJ-1) + X(2,JJ-1)   )
                IF ( MASK(NI,JJ) == 0 )   X(NI,JJ) = &
                     &    0.35*X(NI,JJ) + 0.65*0.125*( &
                     &    X(1+K_EW,JJ)   + X(NI,JJ+1)   + X(NI-1,JJ)   + X(NI,JJ-1) + &
                     &    X(1+K_EW,JJ+1) + X(NI-1,JJ+1) + X(NI-1,JJ-1) + X(1+K_EW,JJ-1)   )
             END IF
             !!
          END DO
          !!
       END DO
       !!
    END IF
    !!
    !!
    DEALLOCATE ( MASKV, MASK2, DATA_OLD )
    !!
    END IF
    !!
END SUBROUTINE DROWN
!!
!!
!!
!!
!!
!!
!!
!!
SUBROUTINE DECISION(IX, IY, NI, NJ, X, MASK, MASK2)
  !!
  !! ARGUMENTS :
  !! -----------
  INTEGER(I4),                   INTENT(IN)    :: IX, IY, NI, NJ
  REAL(R8), DIMENSION(NI,NJ), INTENT(INOUT) :: X
  INTEGER(I4), DIMENSION(NI,NJ), INTENT(IN)    :: MASK
  INTEGER(I4), DIMENSION(NI,NJ), INTENT(INOUT) :: MASK2
  !!
  !!
  !! LOCAL :
  !!--------
  !!  
  INTEGER(I4) :: JI, JJ, JFI, JK
  !!
  REAL(KIND=R8), PARAMETER :: RR=0.707
  REAL(KIND=R8) :: TN, TE, TS, TW, TNE, TSE, TSW, TNW
  INTEGER(I4)      :: WN, WE, WS, WW, WNE, WSE, WSW, WNW
  REAL(KIND=R8) :: RDEN
  !!
  !!
  !!
  !!
  !!     INITIALISATIONS
  !!-----------------------
  JI  = IX
  JFI = IX
  JJ  = IY
  !!
  !!
  !!----------------------------
  !! NORTH :
  !!----------------------------
  IF ( MASK(JI,JJ+1)==1 ) THEN
     TN=X(JI,JJ+1)
     WN=1
  ELSE
     TN=0.
     WN=0
  ENDIF
  !!
  !!----------------------------
  !! NORTH-EAST :
  !!----------------------------
  !! EAST-WEST PERIODICITY :
  IF ( JI == NI ) THEN
     JFI=0
  ELSE
     JFI=JI
  ENDIF
  !!
  IF ( MASK(JFI+1,JJ+1)==1 ) THEN
     TNE=X(JFI+1,JJ+1)
     WNE=1
  ELSE
     TNE=0.
     WNE=0
  ENDIF
  JFI=JI
  !!
  !!----------------------------
  !! EAST :
  !!----------------------------
  !! EAST-WEST PERIODICITY :
  IF ( JI == NI ) THEN
     JFI=0
  ELSE
     JFI=JI
  ENDIF
  !!
  IF ( MASK(JFI+1,JJ)==1 ) THEN
     TE=X(JFI+1,JJ)
     WE=1
  ELSE
     TE=0.
     WE=0
  ENDIF
  JFI=JI
  !!
  !!----------------------------
  ! SOUTH-EAST :
  !!----------------------------
  !! EAST-WEST PERIODICITY :
  IF ( JI == NI ) THEN
     JFI=0
  ELSE
     JFI=JI
  ENDIF
  !!
  IF ( MASK(JFI+1,JJ-1)==1 ) THEN
     TSE=X(JFI+1,JJ-1)
     WSE=1
  ELSE
     TSE=0.
     WSE=0
  ENDIF
  JFI=JI
  !!
  !!----------------------------
  ! SOUTH :
  !!----------------------------
  IF ( MASK(JI,JJ-1)==1 ) THEN
     TS=X(JI,JJ-1)
     WS=1
  ELSE
     TS=0.
     WS=0
  ENDIF
  !!
  !!----------------------------
  ! SOUTH-WEST :
  !!----------------------------
  !! EAST-WEST PERIODICITY :
  IF ( JI == 1 ) THEN
     JFI=NI+1
  ELSE
     JFI=JI
  ENDIF
  !!
  IF ( MASK(JFI-1,JJ-1)==1 ) THEN
     TSW=X(JFI-1,JJ-1)
     WSW=1
  ELSE
     TSW=0.
     WSW=0
  ENDIF
  JFI=JI
  !!
  !!----------------------------
  ! WEST :
  !!----------------------------
  !! EAST-WEST PERIODICITY :
  IF ( JI == 1 ) THEN
     JFI=NI+1
  ELSE
     JFI=JI
  ENDIF
  !!
  IF ( MASK(JFI-1,JJ)==1 ) THEN
     TW=X(JFI-1,JJ)
     WW=1
  ELSE
     TW=0.
     WW=0
  ENDIF
  JFI=JI
  !!
  !!----------------------------
  ! NORTH-WEST :
  !!----------------------------
  !! EAST-WEST PERIODICITY :
  IF ( JI == 1 ) THEN
     JFI=NI+1
  ELSE
     JFI=JI
  ENDIF
  !!
  IF ( MASK(JFI-1,JJ+1)==1 ) THEN
     TNW=X(JFI-1,JJ+1)
     WNW=1
  ELSE
     TNW=0.
     WNW=0
  ENDIF
  JFI=JI
  !!
  !!------------------------
  !!
  !! COMPUTATION OF THE SEA VALUE TO GIVE TO THIS EARTH POINT
  !!----------------------------------------------------------
  RDEN     = (WN+WE+WS+WW+RR*(WNE+WSE+WSW+WNW))
  !!
  X(JI,JJ) = (WN*TN+WE*TE+WS*TS+WW*TW+RR*(WNE*TNE+WSE*TSE+WSW*TSW+WNW*TNW))/RDEN
  !!
  !! FORMER MASK POINT BECOMES SEA POINT :
  MASK2(JI,JJ) = 1
  !!
  !!
END SUBROUTINE DECISION
!!
!!
!!
!!----------------------------------------------------------
!!
!!
!!
!!
SUBROUTINE TEST4SEA(IX, IY, NI, NJ, X, MASK, NDIST, NDIR, NRES)
  !!
  !!#########################################################################
  !!
  !!
  !!  IX  : CURRENT X POSITION                (INTEGER)
  !!  IY  : CURRENT Y POSITION                (INTEGER)
  !!  NI  : X DIMENSION OF ARRAY X            (INTEGER)     
  !!  NJ  : Y DIMENSION OF ARRAY Y            (INTEGER)
  !!  X   : CURRENTLY TREATED ARRAY           (2D ARRAY)
  !!  MASK : MASK ARRAY, 1=SEA, 0=LAND         (2D ARRAY)
  !!  NDIST : DISTANCE (IN POINTS) TO CHECK     (INTEGER)
  !!
  !!                                              1
  !!  NDIR: DIRECTION TO CHECK   (INTEGER)    4 - | - 2
  !!                                              3
  !!
  !!  NRES: RESULT ->  -1 = NO SEA POINT WAS FOUND ON THE WAY NEITHER AT POSITION NDIST 
  !!                  0<N<NDIST = THE FIRST SEA POINT WAS FOUND AT POSITION N
  !!
  !!
  !!#########################################################################
  !!
  INTEGER(I4),                        INTENT(IN) :: IX, IY, NI, NJ, NDIST, NDIR
  INTEGER(I4), DIMENSION(NI,NJ),      INTENT(IN) :: MASK
  REAL(KIND=R8), DIMENSION(NI,NJ), INTENT(IN) :: X
  INTEGER(I4),                        INTENT(OUT):: NRES
  !!
  !!
  !! LOCAL :
  !! -------
  INTEGER(I4) :: JI, JJ, NCPT, JFI
  !!
  JI  = IX
  JFI = IX
  JJ  = IY 
  !!
  !!
  NRES = -1  
  NCPT = 0
  !!
  IF (NDIR==1) THEN
     DO WHILE ((JJ<IY+NDIST).AND.(NRES==-1))
        JJ   = JJ+1
        NCPT = NCPT+1
        IF(MASK(JI,JJ)/=0) THEN
           NRES = NCPT
        ENDIF
     ENDDO
  ENDIF
  !!
  IF (NDIR==2) THEN
     DO WHILE ((JI<IX+NDIST).AND.(NRES==-1))
        IF(JI==NI) JFI=0
        JI   = JI+1
        JFI  = JFI+1
        NCPT = NCPT+1
        IF(MASK(JFI,JJ)/=0) THEN
           NRES = NCPT
        ENDIF
     ENDDO
  ENDIF
  !!
  IF (NDIR==3) THEN
     DO WHILE ((JJ>IY-NDIST).AND.(NRES==-1))
        JJ   = JJ-1
        NCPT = NCPT+1
        IF(MASK(JI,JJ)/=0) THEN
           NRES = NCPT
        ENDIF
     ENDDO
  ENDIF
  !!
  IF (NDIR==4) THEN
     DO WHILE ((JI>IX-NDIST).AND.(NRES==-1))
        IF(JI==1) JFI=NI+1
        JI   = JI-1
        JFI  = JFI-1
        NCPT = NCPT+1
        IF(MASK(JFI,JJ)/=0) THEN
           NRES = NCPT
        ENDIF
     ENDDO
  ENDIF
  !!
  !!
END SUBROUTINE TEST4SEA
!!
!!
!!
!!
!!##########################################################################

 SUBROUTINE FILL_EXTRA_BANDS(K_EW, LX, LY, X, Y, DAT, LXP4, LYP4, XP4, YP4, DATP4)
    !!
    !!============================================================================
    !! EXTENDING INPUT ARRAYS WITH AN EXTRABAND OF TWO POINTS AT NORTH,SOUTH,EAST 
    !! AND WEST BOUNDARIES.
    !!
    !! THE EXTENSION IS DONE THANKS TO AKIMA'S EXPTRAPOLATION METHOD.
    !!
    !! EAST-WEST PERIODICITY OF GLOBAL MAP IS TAKEN INTO ACCOUNT THROUGH 'K_EW' :
    !!
    !!
    !!  K_EW : EAST-WEST PERIODICITY ON THE INPUT FILE/GRID
    !!         K_EW = -1  --> NO PERIODICITY
    !!         K_EW >= 0  --> PERIODICITY WITH OVERLAP OF K_EW POINTS
    !! 
    !!
    !!                       AUTHOR : LAURENT BRODEAU, 2007
    !!============================================================================
    !!
    !!
    INTEGER(I4) ,                       INTENT(IN)  :: K_EW
    !! 
    INTEGER(I4),                        INTENT(IN)  :: LX, LY, LXP4, LYP4
    REAL(R8), DIMENSION(LX,LY),     INTENT(IN)  :: X, Y, DAT
    !!
    REAL(R8), DIMENSION(LXP4,LYP4), INTENT(OUT) :: XP4, YP4, DATP4
    !!
    !! LOCAL
    INTEGER(I4) :: JI, JJ
    !!
    !!
    !!
    !!   C R E A T I N G   E X T E N D E D   A R R A Y S  :
    !!   --------------------------------------------------
    !!
    !! INITIALISING :
    !! --------------
    XP4   = 0.
    YP4   = 0.
    DATP4 = 0.
    !!
    !! FILLING CENTERS :
    !! -----------------
    XP4(3:LXP4-2, 3:LYP4-2)     = X(:,:)
    YP4(3:LXP4-2, 3:LYP4-2)     = Y(:,:)
    DATP4(3:LXP4-2, 3:LYP4-2)   = DAT(:,:)
    !!
    !!
    !! X ARRAY :
    !! ---------
    !!
    IF (K_EW /= -1) THEN   ! WE CAN USE EAST-WEST PERIODICITY OF INPUT FILE TO
       !!                   ! FILL EXTRA BANDS :
       XP4( 1     , 3:LYP4-2) = X(LX - 1 - K_EW , :) - 360.
       XP4( 2     , 3:LYP4-2) = X(LX - K_EW     , :) - 360.
       XP4(LXP4   , 3:LYP4-2) = X( 2 + K_EW     , :) + 360.
       XP4(LXP4-1 , 3:LYP4-2) = X( 1 + K_EW     , :) + 360.
       !!
    ELSE
       !!
       !! LEFT SIDE :
       XP4(2, 3:LYP4-2) = X(2,:) - (X(3,:) - X(1,:))
       XP4(1, 3:LYP4-2) = X(1,:) - (X(3,:) - X(1,:))
       !!
       !! RIGHT SIDE :
       XP4(LXP4-1, 3:LYP4-2) = X(LX-1,:) + X(LX,:) - X(LX-2,:)
       XP4(LXP4  , 3:LYP4-2) = X(LX,:)   + X(LX,:) - X(LX-2,:)
       !!
    END IF
    !!
    !!
    !! BOTTOM SIDE :
    XP4(:, 2) = XP4(:,4) - (XP4(:,5) - XP4(:,3))
    XP4(:, 1) = XP4(:,3) - (XP4(:,5) - XP4(:,3))
    !!
    !! TOP SIDE :
    XP4(:,LYP4-1) = XP4(:,LYP4-3) + XP4(:,LYP4-2) - XP4(:,LYP4-4)
    XP4(:,LYP4)   = XP4(:,LYP4-2) + XP4(:,LYP4-2) - XP4(:,LYP4-4)
    !!
    !!
    !!
    !! Y ARRAY :
    !! ---------
    !!
    !! TOP SIDE :
    YP4(3:LXP4-2, LYP4-1) = Y(:, LY-1) + Y(:,LY) - Y(:,LY-2)
    YP4(3:LXP4-2, LYP4)   = Y(:, LY)   + Y(:,LY) - Y(:,LY-2)
    !! BOTTOM SIDE :
    YP4(3:LXP4-2, 2) = Y(:,2) - (Y(:,3) - Y(:,1))
    YP4(3:LXP4-2, 1) = Y(:,1) - (Y(:,3) - Y(:,1))
    !!
    !!
    IF (K_EW /= -1) THEN   ! WE CAN USE EAST-WEST PERIODICITY
       !!                   ! FILL EXTRA BANDS :
       YP4( 1     , :) = YP4(LX - 1 - K_EW + 2, :)
       YP4( 2     , :) = YP4(LX - K_EW     + 2, :)
       YP4(LXP4   , :) = YP4( 2 + K_EW     + 2, :)
       YP4(LXP4-1 , :) = YP4( 1 + K_EW     + 2, :)
       !!
    ELSE
       !!
       !! LEFT SIDE :
       YP4(2, :) = YP4(4,:) - (YP4(5,:) - YP4(3,:))
       YP4(1, :) = YP4(3,:) - (YP4(5,:) - YP4(3,:))
       !! RIGHT SIDE :
       YP4(LXP4-1,:) = YP4(LXP4-3,:) + YP4(LXP4-2, :) - YP4(LXP4-4, :)
       YP4(LXP4,:)   = YP4(LXP4-2,:) + YP4(LXP4-2,:)  - YP4(LXP4-4, :)
       !!
    END IF
    !!
    !!
    !! DATA ARRAY :
    !! ------------
    !!
    IF (K_EW /= -1) THEN   ! WE CAN USE EAST-WEST PERIODICITY OF INPUT FILE TO
       !!                   ! FILL EXTRA BANDS :
       DATP4( 1     , 3:LYP4-2) = DAT(LX - 1 - K_EW , :)
       DATP4( 2     , 3:LYP4-2) = DAT(LX - K_EW     , :)
       DATP4(LXP4   , 3:LYP4-2) = DAT( 2 + K_EW     , :)
       DATP4(LXP4-1 , 3:LYP4-2) = DAT( 1 + K_EW     , :)
       !!
       !!
    ELSE
       !!
       !! LEFT SIDE :
       DO JJ = 3, LYP4-2
          CALL EXTRA_2_EAST(XP4(LXP4-4,JJ),XP4(LXP4-3,JJ),XP4(LXP4-2,JJ),        &
               &          XP4(LXP4-1,JJ),XP4(LXP4,JJ),                         &
               &          DATP4(LXP4-4,JJ),DATP4(LXP4-3,JJ),DATP4(LXP4-2,JJ),  &
               &          DATP4(LXP4-1,JJ),DATP4(LXP4,JJ) )  
       END DO
       !!
       !! RIGHT SIDE :
       DO JJ = 3, LYP4-2
          CALL EXTRA_2_WEST(XP4(5,JJ),XP4(4,JJ),XP4(3,JJ),                    &
               &          XP4(2,JJ),XP4(1,JJ),                               &
               &          DATP4(5,JJ),DATP4(4,JJ),DATP4(3,JJ),               &
               &          DATP4(2,JJ),DATP4(1,JJ) )  
       END DO
       !!
       !!
    END IF
    !!
    !!
    !! TOP SIDE :
    DO JI = 1, LXP4
       CALL EXTRA_2_EAST(YP4(JI,LYP4-4),YP4(JI,LYP4-3),YP4(JI,LYP4-2),        &
            &          YP4(JI,LYP4-1),YP4(JI,LYP4),                         &
            &          DATP4(JI,LYP4-4),DATP4(JI,LYP4-3),DATP4(JI,LYP4-2),  &
            &          DATP4(JI,LYP4-1),DATP4(JI,LYP4) )  
    END DO
    !!
    !! BOTTOM SIDE :
    DO JI = 1, LXP4
       CALL EXTRA_2_WEST(YP4(JI,5),YP4(JI,4),YP4(JI,3),        &
            &          YP4(JI,2),YP4(JI,1),                    &
            &          DATP4(JI,5),DATP4(JI,4),DATP4(JI,3),    &
            &          DATP4(JI,2),DATP4(JI,1) )  
    END DO
    !!
    !!
  END SUBROUTINE FILL_EXTRA_BANDS
  !!
  !!
  !!
  SUBROUTINE EXTRA_2_EAST(X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5)
    !!
    !!============================================================================
    !!
    !! EXTRAPOLATES 2 EXTRA EAST (OR NORTH) POINTS OF A CURVE WITH AKIMA'S 1D METHOD
    !!
    !! INPUT  : X1, X2, X3, X4, X5, Y1, Y2, Y3
    !! OUTPUT : Y4, Y5
    !!
    !!                       AUTHOR : LAURENT BRODEAU, 2007
    !!============================================================================
    !!
    !!
    REAL(R8), INTENT(IN)  :: X1, X2, X3, X4, X5, Y1, Y2, Y3
    REAL(R8), INTENT(OUT) :: Y4, Y5
    !!
    !! LOCAL :
    REAL(R8) :: A, B, C, D, ALF, BET
    !!
    !!
    A    = X2 - X1
    B    = X3 - X2
    C    = X4 - X3
    D    = X5 - X4
    !!
    ALF  = Y2 - Y1
    BET  = Y3 - Y2
    !!
    IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
       Y4 = Y3 ; Y5 = Y3
    ELSE
       Y4   = C*(2*BET/B - ALF/A) + Y3
       Y5   = Y4 + Y4*D/C + BET*D/B - ALF*D/A - Y3*D/C 
    END IF
    !!
    !!
  END SUBROUTINE EXTRA_2_EAST
  !!
  !!
  !!
  !!
  SUBROUTINE EXTRA_2_WEST(X5, X4, X3, X2, X1, Y5, Y4, Y3, Y2, Y1)
    !!
    !!============================================================================
    !!
    !! EXTRAPOLATES 2 EXTRA WEST (OR SOUTH) POINTS OF A CURVE WITH AKIMA'S 1D METHOD
    !!
    !! INPUT  : X1, X2, X3, X4, X5, Y1, Y2, Y3
    !! OUTPUT : Y4, Y5
    !!
    !!                       AUTHOR : LAURENT BRODEAU, 2007
    !!============================================================================
    !!
    !!
    REAL(R8), INTENT(IN)  :: X1, X2, X3, X4, X5, Y5, Y4, Y3
    REAL(R8), INTENT(OUT) :: Y1, Y2
    REAL(R8) :: A, B, C, D, ALF, BET
    !!
    !! X1 -> X5
    !! X2 -> X4
    !! X3 -> X3
    !! X4 -> X2
    !! X5 -> X1
    !!
    A    = X4 - X5
    B    = X3 - X4
    C    = X2 - X3
    D    = X1 - X2
    !!
    ALF  = Y4 - Y5
    BET  = Y3 - Y4
    !!
    IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
       Y2 = Y3; Y1 = Y3
    ELSE
       Y2   = C*(2*BET/B - ALF/A) + Y3
       Y1   = Y2 + Y2*D/C + BET*D/B - ALF*D/A - Y3*D/C 
    END IF
    !!
    !!
  END SUBROUTINE EXTRA_2_WEST
  !!
  !!
  !!
  !!
END MODULE DRWN
