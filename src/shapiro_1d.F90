      SUBROUTINE SHAPIRO_1D(IM,JM,MASK,RLA_VARIN,ID_NP,RLPA_VAROUT)

      !!=====================================================================
      !!
      !! DESCRIPTION: THIS FUNCTION APPLIES A 1D SHAPIRO FILTER 
      !!              (3 POINTS FILTER) HORIZONTALLY TO A 2D FIELD 
      !!              IN REGULAR GRID
      !! ARGUMENTS :
      !!            RLA_VARIN   : INPUT VARIABLE TO FILTER
      !!            ZLA_MASK    : INPUT MASK VARIABLE
      !!            ID_NP       : NUMBER OF SHAPIRO FILTER ITERATIONS 
      !!            RLPA_VAROUT : OUTPUT FILTERED VARIABLE
      !!
      !! HISTORY : 08/2009  S. CAILLEAU : FROM 1ST VERSION OF N. FERRY
      !!           09/2009  C. REGNIER  : CORRECTIONS
      !!           02/2011  A. STORTO   : ADAPTED FOR OCEANVAR
      !!
      !!=====================================================================
      USE SET_KND
      USE MYFRTPROF
      IMPLICIT NONE
      INTEGER(I4), INTENT(IN)                 :: IM,JM,ID_NP
      REAL(R8), DIMENSION(IM,JM), INTENT(IN)  :: MASK
      REAL(R8), DIMENSION(IM,JM), INTENT(IN)  :: RLA_VARIN
      REAL(R8), DIMENSION(IM,JM), INTENT(OUT) :: RLPA_VAROUT

      REAL(R8), DIMENSION(IM,JM)              :: RLPA_VAROUT_TMP
      REAL(R8), PARAMETER                           :: RL_ALPHA = 0.5_R8
    ! FIXED STABILITY COEFFICIENT (ISOTROPE CASE)
      REAL(R8), PARAMETER                           :: RAP_ANISO_DIFF_XY=2.25_R8
    ! ANISOTROPE CASE
      REAL(R8)                                      :: ALPHAX,ALPHAY, ZNUM, ZDEN,TEST
      INTEGER(I4)                                   :: JI, JJ, JN, NN
!
!! RAP_ANISO_DIFF_XY=2.25 : VALEUR TROUV??E EMPIRIQUEMENT POUR 140 IT??RATION POUR
! LE FILTRE DE SHAPIRO ET 
!! POUR UN RAPPORT D'ANISOTOPIE DE 1.5 : ON FILTRE DE PLUS RAPIDEMENT EN X QU'EN Y.
!
!------------------------------------------------------------------------------
!
! LOOP ON SEVERAL FILTER ITERATIONS
!     GLOBAL OCEAN CASE

             CALL MYFRTPROF_WALL('SHAPIRO_1D: 2D SHAPIRO FILTER',0)

             RLPA_VAROUT(:,:) = RLA_VARIN(:,:)
             RLPA_VAROUT_TMP(:,:) = RLPA_VAROUT(:,:)
!

       ALPHAX=RL_ALPHA
       ALPHAY=RL_ALPHA
!  DX/DY=RAP_ANISO_DIFF_XY  , D_ = VITESSE DE DIFFUSION
!  140 PASSES DU FITRE, LX/LY=1.5, LE RAP_ANISO_DIFF_XY CORRESPONDANT EST:
       IF ( RAP_ANISO_DIFF_XY .GE. 1. ) ALPHAY=ALPHAY/RAP_ANISO_DIFF_XY
       IF ( RAP_ANISO_DIFF_XY .LT. 1. ) ALPHAX=ALPHAX*RAP_ANISO_DIFF_XY

        DO JN = 1,ID_NP   ! NUMBER OF PASSES OF THE FILTER
            DO JI = 2,IM-1
               DO JJ = 2,JM-1
                  ! WE CROP ON THE COAST       
                   ZNUM = RLPA_VAROUT_TMP(JI,JJ)   &
                          + 0.25_R8*ALPHAX*(RLPA_VAROUT_TMP(JI-1,JJ)-RLPA_VAROUT_TMP(JI,JJ))*MASK(JI-1,JJ  )  &
                          + 0.25_R8*ALPHAX*(RLPA_VAROUT_TMP(JI+1,JJ)-RLPA_VAROUT_TMP(JI,JJ))*MASK(JI+1,JJ  )  &
                          + 0.25_R8*ALPHAY*(RLPA_VAROUT_TMP(JI,JJ-1)-RLPA_VAROUT_TMP(JI,JJ))*MASK(JI  ,JJ-1)  &
                          + 0.25_R8*ALPHAY*(RLPA_VAROUT_TMP(JI,JJ+1)-RLPA_VAROUT_TMP(JI,JJ))*MASK(JI  ,JJ+1)
                   RLPA_VAROUT(JI,JJ)=ZNUM*MASK(JI,JJ)+RLA_VARIN(JI,JJ)*(1.-MASK(JI,JJ))
                ENDDO  ! END LOOP JI
            ENDDO  ! END LOOP JJ


!   PERIODICAL CONDITION IN CASE OF CD_OVERLAP (GLOBAL OCEAN)
!   - ON ORCA AND REGULAR GRID WE COPY THE VALUES AT POINTS OF THE PREVIOUS LATITUDE

               RLPA_VAROUT(1,1) = SUM(RLPA_VAROUT(:,2)) / IM
               RLPA_VAROUT(IM,JM) = SUM(RLPA_VAROUT(:,JM-1)) / IM
               RLPA_VAROUT_TMP(:,:) = RLPA_VAROUT(:,:)

         ENDDO  ! END LOOP JN

         CALL MYFRTPROF_WALL('SHAPIRO_1D: 2D SHAPIRO FILTER',1)

END SUBROUTINE SHAPIRO_1D
