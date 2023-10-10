!***  GRID DEFINITIONS
!***

!... Number of Grids supported
INTEGER(I4), PARAMETER :: NGRIDS=7

! Grids resolution in x- and y- direction
REAL(R8), PARAMETER :: ZRESX(NGRIDS) = (/0.5_R8, 0.5_R8, 2._R8, 0.25_R8, 0.0625_R8, 0.037_R8, 1._R8 /)
REAL(R8), PARAMETER :: ZRESY(NGRIDS) = (/0.5_R8, 0.5_R8, 2._R8, 0.25_R8, 0.0625_R8, 0.028_R8, 1._R8 /)

! Start/End of j-index North/South polar regions for ice treatment
INTEGER(I4), PARAMETER,  DIMENSION(NGRIDS) :: NN_END_SOUTHPOLE = (/ 145, 145,  50, 330, 1180, 107, 40 /)
INTEGER(I4), PARAMETER,  DIMENSION(NGRIDS) :: NN_BEG_NORTHPOLE = (/ 310, 310,  100, 670, 2620, 108,141  /)

!... Eof regions mask, if TRUE taken from CNFILEREG
LOGICAL, PARAMETER :: LEOFM(NGRIDS) = (/.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE. /)
LOGICAL, PARAMETER :: LLMDT(NGRIDS) = (/.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE. /)
LOGICAL, PARAMETER :: LLMDTERR(NGRIDS) = (/.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE. /)

!... Grid identification names
CHARACTER(LEN=50), PARAMETER :: CNGRID(NGRIDS) = &
(/ 'GLO05L31 ', 'GLO05L50 ', 'GLO2L31  ', 'GLO025L50', 'GLO116L98', 'BS2736L31', 'WOA1L25' /)

!... File names
CHARACTER(LEN=50), PARAMETER :: CNFILE(NGRIDS) = &
(/ 'GRID.nc', 'GRID.nc', 'GRID.nc', 'GRID.nc', 'GRID.nc', 'GRID_BS.nc', 'GRID.nc' /)

CHARACTER(LEN=50), PARAMETER :: CNFILEREG(NGRIDS) = &
(/ 'GRID_regions.nc', 'GRID_regions.nc', 'GRID_regions.nc', 'GRID_regions.nc', 'GRID_regions.nc', 'GRID_regions.nc', 'GRID_regions.nc' /)

CHARACTER(LEN=50), PARAMETER :: CNFILEMDT(NGRIDS) = &
(/ 'MDT.nc', 'MDT.nc', 'MDT.nc', 'MDT.nc', 'MDT.nc', 'MDT.nc', 'MDT.nc' /)

CHARACTER(LEN=50), PARAMETER :: CNFILEMDTERR(NGRIDS) = &
(/ 'MDT.nc', 'MDT.nc', 'MDT.nc', 'MDT.nc', 'MDT.nc', 'MDT.nc', 'MDT.nc' /)
