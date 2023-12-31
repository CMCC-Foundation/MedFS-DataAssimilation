MODULE OBSDEF

USE SET_KND
IMPLICIT NONE

!... GENERAL
INTEGER(I4), PARAMETER :: NOFAMS=4              ! NUMBER OF OBS FAMILIES
INTEGER(I4), PARAMETER :: NNINS=1              ! INDEX INSITU
INTEGER(I4), PARAMETER :: NNSLA=2              ! INDEX SLA
INTEGER(I4), PARAMETER :: NNSST=3              ! INDEX SST
INTEGER(I4), PARAMETER :: NNSSS=4              ! INDEX SSS

!... DIMENSIONS
INTEGER(I4), PARAMETER :: NOPARAMS=3            ! NUMBER OF PARAMS
INTEGER(I4), PARAMETER :: NOASSIMPARAMS=2       ! NUMBER OF PARAMS
INTEGER(I4), PARAMETER :: NOINSITU=4            ! NUMBER OF INSITU OBS TYPES

!... OBS PARAMETER
INTEGER(I4), PARAMETER :: KKSAL=1               ! SALINITY
INTEGER(I4), PARAMETER :: KKTEMP=2              ! TEMPERATURE
INTEGER(I4), PARAMETER :: KKPOTEMP=3            ! POTENTIAL TEMPERATURE
INTEGER(I4), PARAMETER :: KKSLA=4               ! SEA-LEVEL ANOMALY
INTEGER(I4), PARAMETER :: KKSST=5               ! SEA-SURFACE TEMPERATURE
INTEGER(I4), PARAMETER :: KKSSS=6               ! SEA-SURFACE SALINITY
INTEGER(I4), PARAMETER :: KKPRE=7               ! SEA PRESSURE
INTEGER(I4), PARAMETER :: KKHCDR=8              ! HORIZONTAL CURRENT DIRECTION
INTEGER(I4), PARAMETER :: KKHCSP=9              ! HORIZONTAL CURRENT SPEED
INTEGER(I4), PARAMETER :: KKHCUC=10             ! HORIZONTAL CURRENT U-COMPONENT
INTEGER(I4), PARAMETER :: KKHCVC=11             ! HORIZONTAL CURRENT V-COMPONENT
INTEGER(I4), PARAMETER :: KKSSH=12              ! SEA-SURFACE HEIGHT

!... OBS TYPE
INTEGER(I4), PARAMETER :: KKXBT    = 401        ! XBT/MBT
INTEGER(I4), PARAMETER :: KKTESAC  = 741        ! TESAC
INTEGER(I4), PARAMETER :: KKARGO   = 831        ! ARGO
INTEGER(I4), PARAMETER :: KKBUOY   = 820        ! BUOYS

!... OBS TYPE -- DUMMY VERSION
INTEGER(I4), PARAMETER :: KSXBT    =   1        ! XBT/MBT
INTEGER(I4), PARAMETER :: KSTESAC  =   2        ! TESAC
INTEGER(I4), PARAMETER :: KSARGO   =   3        ! ARGO
INTEGER(I4), PARAMETER :: KSBUOY   =   4        ! BUOYS

CHARACTER(LEN=6), PARAMETER :: CINSNAME(4)=(/'XBT   ','TESAC ','ARGO  ','BUOY  '/)

!***   SATELLITE CONSTANTS

INTEGER(I4), PARAMETER :: NOSLASATS = 9

CHARACTER(LEN=9), PARAMETER :: CSLASATID(NOSLASATS) = (/&
& 'ERS1         ', 'ERS2         ', 'ENVISAT      ',&
& 'GFO          ', 'JASON1       ', 'JASON2       ',&
& 'TP           ', 'CRYOSAT2     ', 'GEOSAT       ' /)

INTEGER(I4), PARAMETER :: CURSAT=999999

INTEGER(I4), PARAMETER :: &
& ISLASATSTART(NOSLASATS)=(/ 199210,199505,200210,200001,200204,200810,199209,201201,198611 /)

INTEGER(I4), PARAMETER :: &
& ISLASATEND  (NOSLASATS)=(/ 199505,200304,CURSAT,200810,CURSAT,CURSAT,200510,CURSAT,198812 /)

INTEGER(I4), PARAMETER :: NOSSTSATS = 6

CHARACTER(LEN=9), PARAMETER :: CSSTSATID(NOSSTSATS) = (/&
& 'AMSR-E       ', 'TMI          ', 'SST MW-OI    ',&
& 'AMSR-E SWATH ', 'TMI SWATH    ', 'REYNOLDS OI  ' /)

INTEGER(I4), PARAMETER :: NOSSSSATS = 2

CHARACTER(LEN=9), PARAMETER :: CSSSSATID(NOSSSSATS) = (/&
& 'AQUARIUS     ', 'SMOS         ' /)

END MODULE OBSDEF
