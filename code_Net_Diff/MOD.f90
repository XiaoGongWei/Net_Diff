! ############Module Setting###############

! =========== Constant variables module =================
module MOD_constant
implicit none
    integer,parameter :: MaxPRN=105  ! max satellite observed at one epoch
    integer :: SatNum=0   ! Satellite Number
    integer :: ParaNum=4   ! Coordinate, Tropsphere parameter, and clock
    integer :: IonoNum=0    ! Ionosphere Number
    integer,parameter :: GNum0=32   ! max GPS satellite number, used for variables defination in MOD
    integer,parameter :: RNum0=27   ! max GLONASS satellite number
    integer,parameter :: CNum0=35   ! max BeiDou satellite number
    integer,parameter :: NumE0=31   ! max GALILEO satellite number
    integer,parameter :: JNum0=7   ! max QZSS satellite number
    integer,parameter :: INum0=7   ! max IRNSS satellite number
    integer :: GNum=32   ! max GPS satellite number
    integer :: RNum=27   ! max GLONASS satellite number
    integer :: CNum=35   ! max BeiDou satellite number
    integer :: NumE=31   ! max GALILEO satellite number
    integer :: JNum=7   ! max QZSS satellite number
    integer :: INum=7   ! max IRNSS satellite number
    real(8), parameter :: c=299792458.0d0   ! velocity of light (m/s)
    ! Elements of WGS84 ellipsolide
    real(8),parameter :: axis=6378137.d0  ! major semi-axis(meter), GGSP standard
    real(8), parameter:: ecc2=0.0066943800229d0  ! square of eccentricity
    real(8), parameter :: pi=3.1415926535897932d0    ! pi
    real(8) :: omg=7.2921151467d-5 !  WGS 84: Earth rotation angular velocity (rad/sec)
    real(8) :: GM=3.986005d14         !  WGS 84: universal gravitational param (m^3/s^2)
                                    !     omg=7.292115d-5      CGCS 2000
                                    !     GM=3.986004418d14  CGCS 2000
    real(8) :: f1     ! frequency
    real(8) :: f2
    real(8) :: f3
    real(8), parameter :: f_L1=10.23e6*154.d0
    real(8), parameter :: f_L2=10.23e6*120.d0
    real(8), parameter :: f_L5=10.23e6*115.d0
    real(8), parameter :: f_B1=10.23e6*152.6d0
    real(8), parameter :: f_B2=10.23e6*118.d0
    real(8), parameter :: f_B3=10.23e6*124.d0
    real(8), parameter :: f_E1=10.23e6*154.6d0
    real(8), parameter :: f_E5a=10.23e6*115.d0
    real(8), parameter :: f_E5b=10.23e6*118.d0
!    f1=10.23d6*152.6d0
!                        f2=10.23d6*118.0d0
!!                        f3=10.23d6*124.0d0
    real(8),parameter :: csGFmax=0.08d0, csGFmin=0.034d0, csGFdt=60.d0
    real(8),parameter :: csMWmax=10.d0, csMWmin=0.9d0, csMWslope=9.d0
    real(8),parameter :: csL1P1maxLength=300.d0,csL1P1max=15.d0, csL1P1init=1.d0, csL1P1slope=5.d0
end module

! =======Set Calculating Time Module=====
module MOD_Time
    implicit none
    type cal_time
       real :: t_begin, t_end
   end type
end module

! =======GLONASS frequency channel=======
module MOD_GLO_Fre
use MOD_constant
    implicit none
    integer(1), save :: Fre_Chann(RNum0)
end module

! =========File Directory Module========
module MOD_FileDir
implicit none
    character(200), save :: ConFile
    character(200), save :: CoorTable
    character(200), save ::  SP3File
    character(200), save ::  ClkFile
    character(200), save ::  EOPFile
    character(200), save ::  AntFile
    character(200), save ::  TideFile  !
    character(200), save ::  PlanetFile
    character(200), save ::  GloFreFile
    character(200), save ::  P1C1File
    character(200), save ::  DCBFile
    character(200), save ::  GPT2GridDir
    character(200), save ::  IONFile
    character(200), save ::  ObsDir
    character(200), save ::  NavDir
    character(200), save ::  sp3Dir
    character(200), save ::  clkDir
    character(200), save ::  OutDir
    character(200), save ::  pcorfile
    character(200), save ::  orbcorrfile
    character(200), save ::  zonecorrfile
    character(200), save ::  BEBFile
    character(200), save ::  IonDir
end module

! =======File ID module=====
module MOD_FileID
implicit none
    integer(1),save :: FileID_Mark=10
    integer(1), save :: ConID
    integer(1), save :: CoorTableID
    integer(1), save :: GLOFreID    ! GLONASS Frequency file(Input)
    integer(1), save :: EOPID    ! Earth Orientation Parameter(Input)
    integer(1), save :: AntID     ! Antenna information, mainly the PCO(Phase Center Offset)(Input)
    integer(2), save :: OLCID   ! Ocean Load Coefficient(Input)
    integer(2), save :: GPTID      ! GPT grid file (input)
    integer(2), save :: P1C1ID      ! P1C1 DCB file (input)
    integer(2), save :: DCBID      ! DCB file (input)
    integer(1), save :: NavID
    integer(1), save :: IONID
    integer(1), save :: IonoID
    integer(1), save :: NavID_R
    integer(1), save :: NavID_C
    integer(1), save :: NavID_E
    integer(1), save :: NavID_J
    integer(1), save :: NavID_I
    integer(1), save :: SP3ID    ! SP3 Data File(Input)
    integer(1), save :: ClkID     ! Clock Data File(Input)
    integer(1), allocatable,save :: ObsID(:)
    integer(1), save :: CoorID
    integer(1), save :: LogID  ! Log File((Output)
     integer(1), save :: CSID   ! For cycle slip 
    integer(1), save :: ResO_CID  ! For Process_Corr
    integer(1), save :: X38ID  ! 
    integer(1), save :: AmbID(91)    ! Ambiguity File(temp)
    integer(1), save :: ModelID(12)   ! used in BRD_SP3
    integer(1), save :: pcorID     ! used in Process_Corr
    integer(1), save :: OrbCorrID   ! used in Process_Corr
    integer(1), save :: BEBID     ! BEB file ID
    integer, save :: LAMBDAID     ! used in double difference
    integer :: amb_success=0, amb_success2=0
end module

! ======EOP module======
module MOD_EOP
implicit none
    type type_EOP
        real(8) :: MJD(2), X(2)=0.d0, Y(2)=0.d0, dUT1(2)=0.d0, dX(2)=0.d0, dY(2)=0.d0
        real(8) :: Xp, Yp
    end type
    type(type_EOP), save :: EOP
    type type_MP
    ! mean pole and their rates at J2000. Values from IERS Conventions (2003),pp84
    ! Used to calculate Pole Tide
    ! Reference :: A_RTK MOD_others.f90
    ! Reference :: http://www.navipedia.net/index.php/Pole_Tide
         real (8) :: xp0 = 0.054d0              ! unit : arcseconds
         real (8) :: yp0 = 0.357d0              ! unit : arcseconds
         real (8) :: xp_rate = 0.00083d0       ! unit : arcseconds per year
         real (8) :: yp_rate = 0.00395d0       ! unit : arcseconds per year
         real (8) :: tref    = 51544.0d0       ! epoch of the reference pole (MJD)
    end type
    type(type_MP), save :: MP
end module

! =====Antenna module======
module MOD_Ant
implicit none
    type type_Ant  ! Phase Center Offset
        character(20) :: Ant_Type
        integer :: PRN=0
        integer :: nzen, nazi
        real(8) :: dazi, dzen
        real(8), allocatable :: Zenith(:)
        real(8), allocatable :: Azimuth(:)
        character(6), allocatable :: Freq(:)
        real(8), allocatable :: PCO(:,:)    ! PCO
        real(8), allocatable :: PCV(:,:,:)
    end type
    type(type_Ant), allocatable ::  Ant(:)
end module

! =====Rotation modile=======
module MOD_Rotation
implicit none
    real(8) :: Rota_S2C(3,3)   ! Rotation from satellite fixed to CRS
    real(8) :: Rota_T2C(3,3)   ! Rotation from TRS to CRS
    real(8) :: Rota_C2T(3,3)   ! Rotation from CRS to TRS
    real(8) :: Rotation(3,3)      ! Appropriate from XYZ to NEU
    real(8)  :: SunCoor(6), MoonCoor(6)
end module

! ======Navigation header module======
module MOD_NavHead
implicit none
    type type_NavHead
        integer :: Version=2
        real(8) :: Alpha(4)=0.d0
        real(8) :: Beta(4)=0.d0
    end type
    type(type_navHead), save :: NavHead
    type(type_navHead), save :: NavHead_R
    type(type_navHead), save :: NavHead_C
    type(type_navHead), save :: NavHead_E
    type(type_navHead), save :: NavHead_J
    type(type_navHead), save :: NavHead_I
end module

module MOD_NavClk
    implicit none
    type type_navclk
        integer :: week(60)
        real(8) :: sow(60)
        real(8) :: a0(60), a1(60)
    end type
    type(type_navclk), save :: NavClk(14)
end module

! ========SP3 File Head module======
module MOD_SP3Head
use MOD_constant
implicit none
    type type_SP3Head
        integer(1) PRNS
        integer(1), allocatable:: PRN(:)
        character(1), allocatable :: SYSTEM(:)
        integer GPSWeek
        real(8) GPSSec
        integer(1) :: OA(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0  ! orbit accuracy
    end type
    type(type_SP3Head), save :: SP3Head
end module

! =========SP3 data module=====
module MOD_SP3Data
use MOD_constant
implicit none
    type type_Eph
        real(8) :: Coor(3,10)=9999.0d0
        real(8) :: Clk(10)=999999.999999d0
    end type
    type type_SP3Data
        integer :: GPSweek(10)=0
        real(8) :: GPSsec(10)=0.d0
        type(type_Eph):: Eph(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)    ! Including GLONASS SP3 Data
    end type
    type(type_SP3Data), save :: SP3Data
end module

! ======Clock Data Module=====
module MOD_ClkData
use MOD_constant
implicit none
    type type_AS
        real(8) :: Clk(2)=9999.d0
        real(8) :: ClkVel=9999.d0
    end type
    type type_ClkData
        integer :: GPSweek(2)
        real(8) :: GPSsec(2)
        type(type_AS) :: AS(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)     ! Including GLONASS Clock data
    end type
    type(type_ClkData), save :: ClkData
end module

! ======Equivalent Satllite Clock Data Module=====
module MOD_ESC
use MOD_constant
implicit none
    real(8) :: ESC(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0.d0
    real(8) :: PCORSow=-1000.d0
    real(8) :: OrbCorr(GNum0+RNum0+CNum0+NumE0+JNum0+INum0,3)=0.d0
    real(8) :: OrbCorrSow=-1000.d0
end module

! =======Navigation data module=======
module MOD_NavData
use MOD_constant
implicit none
    type type_NavPRN
        integer(2) :: GPSweek
        real(8) :: GPSsec
        real(8) :: a0, a1, a2
        real(8)  :: aDot, IODE, Crs, delN, M0         ! line2
        real(8) :: Cuc, e, Cus, sqrtA   ! line3
        real(8) :: toe, Cic, Omega, Cis ! line4
        real(8) :: i0, Crc, w, OmegaDot ! line5
        real(8) :: idot, nDot, Code, WeekNo         ! line 6
        real(8) :: Health, IODC     ! line 7
        real(8) :: TGD(2)  ! TGD1: B1-B3   TGD2: B2-B3
        real(8) :: ISCL1CA=0.d0, ISCL2C=0.d0, ISCL5I5=0.d0, ISCL5Q5=0.d0  ! Line 8, used in CNAV
    end type
    type type_NavData
        type(type_NavPRN),allocatable :: Nav(:)
    end type
    type(type_NavData), save :: NavData(GNum0)
    type(type_NavData), save :: NavData_C(CNum0)
    type(type_NavData), save :: NavData_E(NumE0)
    type(type_NavData), save :: NavData_J(JNum0)
    type(type_NavData), save :: NavData_I(INum0)
end module

! =======GLONASS Navigation data module=======
module MOD_NavData_R
use MOD_constant
implicit none
    type type_NavPRN_R
        integer :: year,mon,day,hour,min
        real(8) :: sec
        integer(2) :: GPSweek
        real(8) :: GPSsec
        real(8) :: a0, a1, a2
        real(8)  :: X,Y,Z
        real(8) :: Vx, Vy, Vz
        real(8) :: Ax, Ay, Az
        real(8) :: Health
    end type
    type type_NavData_R
        type(type_NavPRN_R), allocatable :: Nav(:)
    end type
    type(type_NavData_R), save :: NavData_R(RNum0)
end module

! ===Observation file header module===
module MOD_ObsHead
    type type_ObsHead
        integer(1) :: Version
        integer(1),allocatable :: ObsTypes(:)
        integer(1),allocatable :: ObsType(:,:)
        real(8) :: AppCoor(3)  ! X,Y,Z
        real :: Antenna(3) ! N,E,U
        character(20) :: Ant_Type  ! antenna type
        real(8) :: Interval=30.d0  ! Interval seconds, the default value is 30 s
        integer :: GPSweek  ! Time of first observation
        real(8) :: GPSsec
        integer :: DOY   ! day of year
    end type
    type(type_ObsHead), allocatable, save :: ObsHead(:)
end module

! ======Obs Data module======
module MOD_ObsData
    implicit none
    type type_ObsData
        integer(1) :: Flag
        integer :: year, mon, day, hour, min
        real(8) :: sec
        integer :: week
        real(8) :: sow
        real(8) :: Clk_Bias
        integer(1) :: PRNS
        integer(1), allocatable :: PRN(:)
        character(1), allocatable :: System(:)
        real(8), allocatable :: C1(:), C2(:), P1(:), P2(:),P3(:), L1(:), L2(:), L2C(:), L3(:), D1(:), D2(:), D3(:)
        integer(1), allocatable :: LLI1(:), LLI2(:), LLI2C(:), LLI3(:)   ! Loss of lock indicator, see RINEX Format file
        real(8) :: std(20)
        real(8) :: rela(20)
    end type
end module

!=====Cycle Slip Mark Flag======
module MOD_CycleSlip
use MOD_constant
   implicit none
    type type_CS        
        integer(1) :: arcLengthMW  =  0
        real(8) :: nMWmean=0.d0, nMWmean2=0.d0
        integer(1) :: arcLengthGF  =  0
        integer :: WeekPrev(3) = 0
        real(8) :: SowPrev(3) = 0.d0, GFPrev(3)
        integer(2) :: arcLengthL1P1  =  0
        integer :: LastWeek
        real(8) :: LastSow
        real(8) :: L1P1mean, L1P1mean2
        real(8) :: OMC(4)=0.d0
        integer(1)  :: Slip        =  0
    end type
   type type_CycleSlip
!       real(8) :: GF(2,GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=9999.0d0
!       real(8) :: MW(2,GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=9999.0d0
!       real(8) :: dT(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0.d0
!       integer :: n(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0
!       real(8) :: P1(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0.d0, P2(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0.d0
!       real(8) :: PreL1(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0.d0, PreL2(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0.d0
       real(8) :: OMC(GNum0+RNum0+CNum0+NumE0+JNum0+INum0,4) = 0.d0  ! For satellite difference cyccle slip detection
       integer(1)  :: Slip(GNum0+RNum0+CNum0+NumE0+JNum0+INum0) =  0  ! For satellite difference cyccle slip detection
       integer(1) :: CScount=99  ! For satellite difference cyccle slip detection
       type(type_CS) :: CS(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)   ! Cycle slip
   end type
   type(type_CycleSlip),allocatable ::  CycleSlip(:)
end module 

! ========= Station module ==========
module MOD_STA
use MOD_constant
    type type_Trop
        real(8) :: press,temp, dT, Tm, e, ah, aw, la, undu
        real(8) :: map_dry, map_wet
        real(8) :: ZHD, ZWD
    end type
    type type_Pre
        real(8) :: n=0.d0
        real(8) :: CBias1, CBias2, Mp1=0.d0, Mp2=0.d0
        real(8) :: Range=0.d0, PrePhase=0.d0, HisPhase(60)=99.d0,HisRange(60)=99.d0
        real(8) :: n1=0.d0, n2=0.d0, P1=0.d0, P2=0.d0, DP=0.d0, Sow1=0.d0, Sow2=0.d0
    end type
    type type_OLC
        real(8) :: OLC(11,6)
        real(8) :: CMC(11,6)
        logical :: found
    end type
    type type_STA_STA
        character(4) :: Name ! Station Name
        character(4) :: SKD  ! Station coordinate status, fix or not
        real(8) :: Coor(3) ! Approximate Station Coordinate
        real(8) :: TrueCoor(3) ! True Station Coordinate
        real(8) :: XYZ(3) ! Positioning coordiante of each epoch
        logical :: flag_InitialCoor=.false. ! Flag of initial Station Coordinate
        real(8) :: BLH(3) ! Station BLH
        real(8) :: NEU(3) ! Antenna
        character(20) :: Ant_Type  ! antenna type
        real(8) :: Rotation(3,3)
        type(type_Trop) :: Trop
        type(type_Pre) :: Pre(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)   ! Pseudorange smooth
        type(type_OLC) :: OLC    ! ocean load coefficient
    end type
    type type_STA
        integer :: Num=0 ! total station number
        integer :: FixNum=0   ! number of fixed stations 
        type(type_STA_STA), allocatable :: STA(:)
    end type
    type(type_STA), save :: STA
end module

! ======= Residual Information module ========
! used in station difference
module MOD_Res
use MOD_constant
    type type_Res
        integer :: N=0
        real(8) :: sow=-1000.d0
        integer :: PRN(MaxPRN*2)=0
        character(2) :: Code(MaxPRN*2)=""
        real(8) :: L(MaxPRN*2)=0.d0
        integer :: CS(MaxPRN*2)=0
    end type
     type(type_Res),  allocatable, save :: Res(:)
end module

! ========== Zero Difference module =========
! used in double difference
module MOD_ZD
use MOD_constant
! integer, parameter :: MaxPRN=35
    type type_ZD
        integer :: week                 =  0
        real(8)  :: sow                   =  0.d0
        integer :: PRNS                 =  0
        integer(1) :: Sys(MAXPRN)=0
        character(1) :: System(MAXPRN)
        integer :: PRN(MAXPRN)  =  0
        integer(1) :: PRN_S(MAXPRN) =  0
        real(8)  :: Ele(MAXPRN)    =  0.d0
        real(8)  :: P(MAXPRN)       =  0.d0
        real(8), allocatable  :: A(:, :)
        real(8)  :: Corr(MAXPRN)  =  0.d0
        real(8)  :: s(MAXPRN)       =  0.d0
        real(8)  :: P1(MaxPRN)     =  0.d0
        real(8)  :: P2(MaxPRN)     =  0.d0
        real(8)  :: P1CS(MaxPRN)     =  0.d0   ! for cycle slip
        real(8)  :: P2CS(MaxPRN)     =  0.d0
        real(8)  :: amb0(MaxPRN,3)     =  0.d0
        real(8)  :: L1(MaxPRN)     =  0.d0
        real(8)  :: L2(MaxPRN)     =  0.d0
        real(8)  :: WL(MaxPRN)    =  0.d0
        real(8)  :: W4(MaxPRN)    =  0.d0
        real(8)  :: EWL(MaxPRN)    =  0.d0
        real(8)  :: WL_amb(MaxPRN)    =  0.d0
        real(8)  :: WL_amb_n(MaxPRN)    =  0.d0
        real(8)  :: EWL_amb(MaxPRN)    =  0.d0
    end type
!    type(type_ZD) :: ZD(2)  ! 1: Reference station; 2: User station
end module

! ========== Station Difference module =========
! used in double difference
module MOD_SD
use MOD_constant
! integer, parameter :: MaxPRN=35
    type type_SD
        integer :: week                 =  0
        real(8)  :: sow                   =  0.d0
        integer :: PRNS                 =  0
        integer(1) :: Sys(MAXPRN)=0
        character(1) :: System(MAXPRN)
        integer :: PRN(MAXPRN)  =  0
        integer(1) :: PRN_S(MAXPRN) =  0
        real(8)  :: Ele(MAXPRN)    =  0.d0
        real(8)  :: Q(MAXPRN)      =  0.d0
        real(8), allocatable  :: A(:, :)
        real(8)  :: Corr(MAXPRN)  =  0.d0
        real(8)  :: s(MAXPRN)       =  0.d0
        real(8)  :: P1(MaxPRN)     =  0.d0
        real(8)  :: P2(MaxPRN)     =  0.d0
        real(8)  :: P1CS(MaxPRN)     =  0.d0   ! for cycle slip
        real(8)  :: P2CS(MaxPRN)     =  0.d0
        real(8)  :: L1(MaxPRN)     =  0.d0
        real(8)  :: L2(MaxPRN)     =  0.d0
        real(8)  :: WL(MaxPRN)    =  0.d0
        real(8)  :: W4(MaxPRN)    =  0.d0
        real(8)  :: WL_amb(MaxPRN)    =  0.d0
        real(8)  :: EWL(MaxPRN)    =  0.d0
        real(8)  :: EWL_amb(MaxPRN)    =  0.d0
    end type
!    type(type_SD), save :: SD
end module

! ========== Double Difference module =========
! used in double difference
module MOD_DD
use MOD_constant
! integer, parameter :: MaxPRN=35
    type type_DD
        integer :: week                =  0
        real(8)  :: sow                  =  0.d0
        integer :: PRNS                =  0     ! not include the reference satellite
        integer(1) :: Sys(MAXPRN)=0
        character(1) :: System(MAXPRN)
        integer :: RefSat(5)          =  0     ! reference satellite
        integer :: PRN(MAXPRN) =  0
        integer(1) :: PRN_S(MAXPRN) =  0
        real(8)  :: Ele(MAXPRN)    =  0.d0
        real(8)  :: P(MaXPRN, MaxPRN)       =  0.d0
        real(8)  :: Q(MaxPRN,MaxPRN)       =  0.d0
        real(8), allocatable  :: A(:, :)
        real(8)  :: Corr(MAXPRN)  =  0.d0
        real(8)  :: s(MAXPRN)       =  0.d0
        real(8)  :: P1(MaxPRN)     =  0.d0
        real(8)  :: P2(MaxPRN)     =  0.d0
        real(8)  :: L1(MaxPRN)     =  0.d0
        real(8)  :: L2(MaxPRN)     =  0.d0
        real(8)  :: WL(MaxPRN)    =  0.d0
        real(8)  :: W4(MaxPRN)    =  0.d0
        real(8)  :: WL_amb(MaxPRN)    =  0.d0
        real(8)  :: EWL(MaxPRN)    =  0.d0
        real(8)  :: EWL_amb(MaxPRN)    =  0.d0
    end type
!    type(type_DD), save :: DD
end module

! ========== Previous Station Difference Information ==========
! used in cycle slip detection
module MOD_PreSD
use MOD_constant
! integer, parameter :: MaxPRN=35
    type type_PreSD
        integer :: WeekGF      =   0
        real(8)  :: sowGF         =   0.d0
        real(8) :: GF                =   0.d0
        integer :: WeekMW   =   0
        real(8)  :: sowMW      =   0.d0
        real(8) :: MW             =   0.d0
        integer(1) :: LLI1=0, LLI2=0, LLI3=0

        integer(1) :: arcLengthMW  =  0
        real(8) :: nMWmean=0.d0, nMWmean2=0.d0
        integer(1) :: arcLengthGF  =  0
        integer :: WeekPrev(3) = 0
        real(8) :: SowPrev(3) = 0.d0, GFPrev(3)
        integer(2) :: arcLengthL1P1  =  0
        integer :: LastWeek
        real(8) :: LastSow
        real(8) :: L1P1mean, L1P1mean2
    end type
    type(type_PreSD) :: PreSD(MaxPRN)
end module

! ========== Normal equation module =========
! used in NEQ solution
module MOD_NEQ
use MOD_constant
! integer, parameter :: MaxPRN=35
    type type_NEQ
        integer :: PRNS=0   ! PRNS in this epoch
        integer  :: Npar =  0
        integer :: SumN     ! total error equations
!        real(8)  :: A((MaxPRN-10)*4,MaxPRN*2+3), L((MaxPRN-10)*4), V((MaxPRN-10)*4)
        integer :: N=MaxPRN*2+3  ! parameter numbers
        integer  :: PRN(MAXPRN) =  0
        integer(1) :: Sys(MAXPRN)=0
        character(1) :: System(MAXPRN)
        real(8)  :: Ele(MAXPRN)    =  0.d0
        real(8)  :: R(MaxPRN, MaxPRN)=0.d0
        real(8)  :: P(MaxPRN, MaxPRN)=0.d0
        real(8), allocatable :: Ap1(:, :), Ap2(:, :), Awl(:, :), Aw4(:, :), Aewl(:, :)
        real(8)  :: Lp1(MaxPRN), Vp1(MaxPRN)
        real(8)  :: Lp2(MaxPRN), Vp2(MaxPRN)
        real(8)  :: Lwl(MaxPRN), Vwl(MaxPRN)
        real(8)  :: Lw4(MaxPRN), Vw4(MaxPRN)
        real(8)  :: Lewl(MaxPRN), Vewl(MaxPRN)
        real(8)  :: maxV(5)   ! 记录每种组合误差方程最大的的残差
        integer :: maxL(5)    ! 记录每种组合最大的残差所在的位置
        real(8), allocatable  :: Nbb(:, :)
        real(8), allocatable  :: InvN(:, :)
        real(8), allocatable  :: U(:)
        real(8), allocatable  :: dx(:)
        real(8) :: amb_WL(MaxPRN)=0.d0
        real(8) :: amb_W4(MaxPRN)=0.d0
        real(8)  :: amb_EWL(MaxPRN) = 0.d0
        real(8) :: fixed_amb(MaxPRN*2)=0.99d0
        integer :: fixed_amb_num(MaxPRN*2)=0
        real(8) :: fixed_amb_ele(MaxPRN*2)=0.d0
        real(8) :: iono(MaxPRN)=0.d0
        real(8) :: amb_L1(MaxPRN)=0.d0
        real(8) :: amb_L2(MaxPRN)=0.d0
        real(8) :: ratio=0.d0
    end type
end module

! ========== Normal equation module =========
module MOD_Epo_NEQ
use MOD_constant
implicit none
    type type_Epo_NEQ
        integer :: PRNS=0
        integer :: SumN
!        real(8)  :: A((MaxPRN-10)*4,3), L((MaxPRN-10)*4), V((MaxPRN-10)*4)
        integer :: N=3  ! parameter numbers
        integer  :: PRN(MAXPRN) =  0
        integer(1) :: Sys(MAXPRN)=0
        character(1) :: System(MAXPRN)
        real(8)  :: Ele(MAXPRN)    =  0.d0
        real(8)  :: R(MaxPRN, MaxPRN)=0.d0
        real(8)  :: P(MaxPRN, MaxPRN)=0.d0
        real(8), allocatable :: Ap1(:, :),Ap2(:, :), Al1(:, :), Al2(:, :), Awl(:, :), Aw4(:, :)
        real(8)  :: Lp1(MaxPRN), Vp1(MaxPRN)
        real(8)  :: Lp2(MaxPRN), Vp2(MaxPRN)
        real(8)  :: Ll1(MaxPRN), Vl1(MaxPRN)
        real(8)  :: Ll2(MaxPRN), Vl2(MaxPRN)
        real(8)  :: Lwl(MaxPRN), Vwl(MaxPRN)
        real(8)  :: Lw4(MaxPRN), Vw4(MaxPRN)
        integer :: ParPRN=0 ! Partial AR PRN
        real(8)  :: maxV(6)=0.d0   ! 记录每种组合误差方程最大的的残差
        integer :: maxL(6)=0.d0    ! 记录每种组合最大的残差所在的位置
        real(8), allocatable  :: Nbb(:,:)   ! for partial ambiguity resolution, only one satellite currently
        real(8), allocatable  :: InvN(:,:)
        real(8), allocatable  :: U(:)
        real(8), allocatable :: dx(:)
        real(8) :: amb_L1(MaxPRN)=0.d0
        real(8) :: amb_L2(MaxPRN)=0.d0
        real(8) :: amb_WL(MaxPRN)=0.d0  ! This is only for wide lane ambiguity rounding. Just for test.
        real(8) :: amb_W4(MaxPRN)=0.d0
        real(8) :: fixed_amb(MaxPRN*2)=0.99d0
        integer :: fixed_amb_num(MaxPRN*2)=0
        real(8) :: ratio=0.d0
    end type
end module

! ========== Doppler Normal Equation module ============
module MOD_NEQ_DP
implicit none
    type type_NEQ_DP
        integer :: PRNS
        integer :: PRN(40)
        real(8) :: Ele(40), P(40)
        real(8) :: A(40,4), L(40)
        real(8) :: Nbb(4,4), InvN(4,4), U(4)
        real(8) :: dx(4), V(40)
        real(8) :: sigma0
        real(8)  :: maxV   ! 记录最大的的残差
        integer :: maxL   ! 记录最大的残差所在的位置
        real(8) :: Sow(5), dt
        real(8) :: Vel(3,5)
        integer(1) :: Flag_Sln(5)
        real(8) :: Coor(3)
    end type
    type(type_NEQ_DP) :: NEQ_DP
end module

! ========== Ionosphere residual module =========
module MOD_IonoDDRes
use MOD_constant
implicit none
    type type_IonoDDRes
        real(8) :: Sow(5)
        real(8) :: L1(5)=0.d0
        real(8) :: L2(5)=0.d0
        real(8) :: amb1(5)=0.d0
        real(8) :: amb2(5)=0.d0
        real(8) :: dL1=0.d0
        real(8) :: dL2=0.d0
    end type
    type (type_IonoDDRes) :: IonoDDRes(MAXPRN)
    integer :: Last_RefSat(5) =  0     ! reference satellite
end module

  ! =======Variable module======
 module MOD_Var
use MOD_constant
 implicit none
     integer :: int_year, int_doy
    character(8) :: str_ymd
     real(8) :: Interval=30.d0
     real(8) :: LimEle=10.d0
     real(8) :: LimSNR=0.d0
     real(8) :: MaxPDOP=100.d0
     real(8) :: Arclength=0.d0
     integer :: arc_epochs
     real(8) :: TropLen=7200.d0
     real(8) :: delay0=0.d0
     real(8) :: delay=0.d0
     character(1) :: Pos_State="S"
     character(2) :: diffmod="DD"
     integer(1)  :: proc_mod=0
     character(4) :: freq_comb="L1L2"
     logical :: Combination(3)=.false.   ! PC/LC/Doppler
     real(8)  ::  sigPC=1.d0,sigLC=0.01d0, sigDP=0.1d0
     logical :: SystemUsed(6)=.false.   ! G/R/C/E/J/I
     integer(1) :: INT_SystemUsed(6)=0   ! G/R/C/E/J/I
     character(2) :: ADmethod="LS"
     logical :: If_ISB=.true.
     logical :: If_IFB=.false.
     character(2) :: ISB_Mode='RW'
     character(2) :: IFB_Mode='FD'
     real(8) :: ISB_step=1.d0
     integer :: GPSweek_st=0, GPSweek_end=3000
     real(8) :: GPSsec_st=0.d0, GPSsec_end=0.d0
     character(2) :: ObsCombine="PC"
     logical :: If_Est_Iono=.false.
     character(1) :: Var_smooth="n"
     character(7) :: cdattype='GPT'    ! meteorological parameter type
     character(5) :: cztd='SAAS'     ! Zenith dry and wet delay model
     character(5) :: cmap='NMF'    ! Mapping function
     character(20) :: CSmethod="GF+MW"   ! cycle slip detect method
     character(5) :: Smooth_Method
     character(4) :: Smooth_Combine='P1P2'
     real(8) :: Smooth_Time
     character(3) :: ObsType='RNX'
     integer(1) :: clktype=1
     integer(1) :: IorQ=0   ! only for X71
     character(3) :: Orbit
     character(5) :: Clk
     logical(1) :: IsCNAV=.false.
     integer(1) :: SatSelected(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=1
     real(8) :: Leap_sec=0.d0
     character(3) :: TimeSys=""
     real(8) :: ObsTime=0.d0,SP3Time=0.d0
     real(8) :: DCB(GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0.d0
     real(8) :: DCBBSX(5,GNum0+RNum0+CNum0+NumE0+JNum0+INum0)=0.d0
     real(8) :: DCBIQ(5,35)
     real(8) :: BEB(3,35)=0.d0  ! BeiDou Broadcast Ephemeris Bias
     real(8) :: iono_14para_day(13,14)
     integer(1) ::  iontype
    ! The following parameter is used in double difference
    integer(1) :: ar_mode=1
    real(8) :: FixEle=20.d0, HoldEle=30.d0  ! Fix and hold satellite elevation
    real(8) :: minratio=3.d0
    logical(1) :: partial_AR=.true.
    integer(1) :: parARnum=2
    logical :: If_Est_WL=.false.   ! This is only for long baseline RTK
    real(8) :: Baseline, Diff_Hgt, Min_Lat   ! This is only for long baseline RTK
    integer :: GloParaNum=0   ! This is for GLONASS float RTK
    real(8) ::   a1= 1.d0 , a2= 0.d0 ! L1
    real(8) ::  b1= 0.d0 , b2= 1.d0  ! L2
    logical :: If_IonoCompensate=.false.
     integer(1) :: Vel_Used=0  ! 0: no velocity; 1: Doppler; 2: IMU
!    real(8),parameter ::   a1= 1.d0 , a2= 0.d0 ! L1 
!    real(8), parameter ::  b1= 0.d0 , b2= 1.d0  ! L2 
!    real(8),parameter ::   a1= 1.d0 , a2= -1.d0 ! Wide-Lane
!    real(8), parameter ::   b1=4.d0, b2= -5.d0  ! W4
!    real(8), parameter ::   b1= -1.d0, b2=2.d0  ! W1
!    real(8), parameter ::   a1= -2.d0, a2=3.d0  ! W2
!    real(8), parameter ::   b1= -3.d0, b2=4.d0  ! W3
!    real(8), parameter ::   a1=154.d0/34.d0, a2= -120.d0/34.d0 ! Iono-free
!    real(8), parameter ::   a1=152.6d0/(152.6d0-118.0d0), a2= -118.d0/(152.6d0-118.0d0) ! Iono-free B1B2
!    real(8), parameter ::   a1=152.6d0/(152.6d0-124.0d0), a2= -124.d0/(152.6d0-124.0d0) ! Iono-free B1B3
!    real(8),parameter ::   b1= 0.d0 , b2= 0.d0 ! Iono-free
 end module