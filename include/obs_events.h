!... Events -- Observations pre-processing
!    A.Storto
!    J.Pistoia ad too shallow obs

INTEGER(i4), PARAMETER :: knoeve=32       ! Number of Events accounted for

INTEGER(i4), PARAMETER :: keve_RETA=0     ! Retained Observation
INTEGER(i4), PARAMETER :: keve_INTE=1     ! Interpolation problems
INTEGER(i4), PARAMETER :: keve_TIME=2     ! Out of assimilation time-window
INTEGER(i4), PARAMETER :: keve_MDTU=3     ! MDT not available (SLA only)
INTEGER(i4), PARAMETER :: keve_BGQC=4     ! Background Quality Check
INTEGER(i4), PARAMETER :: keve_LONM=5     ! NoMotion Lev unapplicable (SLA only)
INTEGER(i4), PARAMETER :: keve_GEOR=6     ! Geographical rejection
INTEGER(i4), PARAMETER :: keve_HTH1=7     ! Horiz. Thinning rejection
INTEGER(i4), PARAMETER :: keve_HTH2=8     ! Horiz. Thinning rejection
INTEGER(i4), PARAMETER :: keve_VTHN=9     ! Vert. Thinning rejection
INTEGER(i4), PARAMETER :: keve_BLAC=10    ! Blacklisted Source
INTEGER(i4), PARAMETER :: keve_DEPL=11    ! Unsupported depth
INTEGER(i4), PARAMETER :: keve_MASK=12    ! Mask problems
INTEGER(i4), PARAMETER :: keve_TOOD=13    ! Too deep observation
INTEGER(i4), PARAMETER :: keve_DSTC=14    ! Obs too close to coast
INTEGER(i4), PARAMETER :: keve_BGCH=15    ! Background value not valid
INTEGER(i4), PARAMETER :: keve_ICER=16    ! Obs below model sea-ice
INTEGER(i4), PARAMETER :: keve_ICED=17    ! Obs too close to sea-ice
INTEGER(i4), PARAMETER :: keve_AGGT=18    ! In-situ aggressive thinning
INTEGER(i4), PARAMETER :: keve_AOER=19    ! Obs error failed
INTEGER(i4), PARAMETER :: keve_TDST=20    ! Obs out of time-window
INTEGER(i4), PARAMETER :: keve_VCHK=21    ! Insitu Obs with rejected obs above
INTEGER(i4), PARAMETER :: keve_CLIM=22    ! Too large departure from climatology
INTEGER(i4), PARAMETER :: keve_REDC=23    ! Redundancy check (Argo vs GTSPP)
INTEGER(i4), PARAMETER :: keve_PDOM=24    ! Partial Domain
INTEGER(i4), PARAMETER :: keve_SBCR=25    ! Bias Correction Failed
INTEGER(i4), PARAMETER :: keve_NSST=26    ! Partial Rejection in case of 2-step SST assimilation
INTEGER(i4), PARAMETER :: keve_UDEN=27    ! Partial Rejection in case of 2-step SST assimilation
INTEGER(i4), PARAMETER :: keve_SSSR=28    ! SSS Specific
INTEGER(i4), PARAMETER :: keve_UNSP=29    ! Parameter Not assimilated
INTEGER(i4), PARAMETER :: keve_REGR=30    ! Regional rejection
INTEGER(i4), PARAMETER :: keve_TOOS=31    ! Too SHALLOW observation
INTEGER(i4), PARAMETER :: keve_TFEW=32    ! Too few observations

INTEGER(i4), PARAMETER :: max_subgroup = 10
INTEGER(i4),DIMENSION(0:knoeve,max_subgroup) :: kcount=0
INTEGER(i4),DIMENSION(0:knoeve) :: kcountt=0

CHARACTER(LEN=42) :: CSTREVE(0:knoeve)

CSTREVE(0) = 'Retained Observation                      '
CSTREVE(1) = 'Interpolation problems                    '
CSTREVE(2) = 'Out of assimilation time-window           '
CSTREVE(3) = 'MDT not available (SLA only)              '
CSTREVE(4) = 'Background Quality Check                  '
CSTREVE(5) = 'Level-of-no-Motion unapplicable (SLA only)'
CSTREVE(6) = 'Geographical rejection                    '
CSTREVE(7) = 'Horizontal thinning rejection, step 1     '
CSTREVE(8) = 'Horizontal thinning rejection, step 2     '
CSTREVE(9) = 'Vertical thinning rejection (IN-SITU only)'
CSTREVE(10)= 'Blacklisted Source                        '
CSTREVE(11)= 'Unsupported depth                         '
CSTREVE(12)= 'Mask problems                             '
CSTREVE(13)= 'Too deep observation                      '
CSTREVE(14)= 'Observation too close to coast            '
CSTREVE(15)= 'Background value not valid                '
CSTREVE(16)= 'Observation below model sea-ice           '
CSTREVE(17)= 'Observation too close to sea-ice          '
CSTREVE(18)= 'Aggressive thinning (IN-SITU only)        '
CSTREVE(19)= 'Failure of obs error specification        '
CSTREVE(20)= 'Out of assimilation time-window           '
CSTREVE(21)= 'In-situ Observation with rejections above '
CSTREVE(22)= 'Too large departure from climatology      '
CSTREVE(23)= 'Redundancy check (IN-SITU only)           '
CSTREVE(24)= 'Partial domain check                      '
CSTREVE(25)= 'Bias Correction Failed                    '
CSTREVE(26)= 'Rejection in case of 2-step Assimilation  '
CSTREVE(27)= 'Undefined background density level        '
CSTREVE(28)= 'SSS Specific Rejections                   '
CSTREVE(29)= 'Parameter rejected for assimilation       '
CSTREVE(30)= 'Regional rejection (ad hoc studies)       '
CSTREVE(31)= 'Too Shallow observation                   '
CSTREVE(32)= 'Too few observations                      '

