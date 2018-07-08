! ======================  Cal_Sat_Pos  ====================
!
! PURPOSE:
!             Calculate the Satellite Coordinate, relativity error, 
!    the range between the satellite and the receiver, using the precise orbit
!
! REFERENCES:
!     For satellite coordinates:
!              http://www.navipedia.net/index.php/Emission_Time_Computation
!               http://www.navipedia.net/index.php/Precise_GNSS_Satellite_Coordinates_Computation
!     For relativity:
!               http://www.navipedia.net/index.php/Relativistic_Clock_Correction
!               http://www.navipedia.net/index.php/Relativistic_Path_Range_Effect
!              http://www.extinctionshift.com/SignificantFindings06B.htm
!
! INPUTS:
!    GPSweek           GPS week
!    GPSsec              GPS second
!    PRN                   satellite number, integer
!    Rec_Coor           coordinate of the station, unit in meter
!    t1                       signal transfer time using LC combined range, unit in second
!    Rela_flag            relativity flag (logical)
!                              .true. : calculate relativity
!                              .false. :  relativity=0
!
! OUTPUTS:
!    Sat_Coor2          satellite coordinate, unit in meter
!    R                        range between the satellite and the receiver, unit in meter
!    Rela                   relativity correction, unit in meter
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================End of Header==============

subroutine Cal_Sat_Pos(GPSweek, GPSsec, PRN, Rec_Coor,t1, Rela_Flag,Sat_Coor2, Sat_Vel, R, Rela, Sat_Clk)
use MOD_constant
use MOD_SP3Data
implicit none
    ! Intent in
    integer ::  PRN
    integer :: GPSweek
    real(8) :: GPSsec
    real(8) :: Rec_Coor(3), t1
    logical :: Rela_Flag
    ! Intent out
    real(8) ::  Sat_Coor(3), R, Sat_Coor2(3)
    real(8) ::  Rela, Sat_Clk
    
    real(8) :: TT(10), t2
    real(8) :: Rotation(3,3)
    real(8) :: TT2(10), VCoor(3,9), Sat_Vel(3)
    real(8) :: Clk(2), TTT(2)
    integer(1) :: i, j
    integer::k
    real(8) :: r_sta, r_sat, r_sta2sat, dbending
    
    TT=(SP3Data%GPSweek-GPSweek)*604800.0d0+SP3Data%GPSsec-GPSsec
    if ( all(dabs(SP3Data%Eph(PRN)%Coor-9999.0d0)<0.1d0) ) then
        Sat_Coor2=9999.0d0
        return
    end if
    
    do while(.true.)
        call Lagrange(TT,SP3Data%Eph(PRN)%Coor,-t1,Sat_Coor,10)  ! Lagrange Interpolation (10 points)
        
        Sat_Coor=Sat_Coor*1000.0d0   ! Unit in meter
        ! Earth rotation correction
        Rotation=reshape( (/ dcos(omg*t1), dsin(omg*t1), 0.d0, -dsin(omg*t1),dcos(omg*t1), 0.d0, 0.d0, 0.d0, 1.0d0 /), (/3,3/) )
        !Rotation= reshape( (/ 1D0, dsin(omg*t1), 0.d0, -dsin(omg*t1),1.d0, 0.d0, 0.d0, 0.d0, 1.d0 /), (/3,3/) )
        
        Sat_Coor2=MATMUL(Sat_Coor,Rotation)
        ! Diatance from the satellite to the reciver
        R=dsqrt(DOT_PRODUCT((Sat_Coor2-Rec_Coor),(Sat_Coor2-Rec_Coor)))
        ! Signal transfering time
        t2=R/c
        if (dabs(t2-t1)<1.0D-10) exit
        t1=t2
    end do


    Rela=0.d0
    if (Rela_Flag) then
        call Lagrange_Vel(TT,SP3Data%Eph(PRN)%Coor,-t1,Sat_Vel,10)
        ! **** Relativistic Clock Correction *****
        !    The clocks, one placed in the satellite and the other ion the terrestrial surface, will differ due to 
        ! the difference of gravitational potential and speed between them.
        !    A constant componant has already modified (in factory) the clockoscillating frequency of the satellite.
        !    A periodical componant due to orbit eccentricity must be corrected by the user receiver software.
        Sat_Vel=Sat_Vel*1000.d0
        Rela=-2.0d0/c*(DOT_PRODUCT(Sat_Coor,Sat_Vel)) ! Unit in meter ! 见 GPS ICD
        Sat_Vel=(Sat_Coor-Sat_Coor2)/t1+Sat_Vel
        ! **** Relativistic Path Range Effect *****
        !     This is a secondary relativistic effect that can be required only for high accuracy positioning.
        ! Its net effect on range is less than 2 cm and thence, for most purpose it can be neglected.
        !      This effect is named the Shapiro signal propagation delay and introduces a general relativistic 
        !  correction to the geometric range. The Shapiro delay was first noticed by Irvin L. Shapiro in 1964,
        !  the transit time required for a microwave signal to propagate through space, which was affected
        !  by the relative position of the sun.( http://www.extinctionshift.com/SignificantFindings06B.htm)
        !     In this programe, the sign of dbending is oppsite with relativistic clock correction.
        ! 实际上，Relativistic Path Range Effect变化不大，RMS只有2mm左右，对定位影响只有4mm左右。
        if (all(Rec_Coor/=0.d0)) then
            r_sta=dsqrt( DOT_PRODUCT(Rec_Coor,Rec_Coor) )
            r_sat=dsqrt( DOT_PRODUCT(Sat_Coor,Sat_Coor) )
            r_sta2sat=dsqrt( DOT_PRODUCT( (Sat_Coor-Rec_Coor),(Sat_Coor-Rec_Coor) ) )
            dbending=2.d0*GM/c**2*log((r_sta+r_sat+r_sta2sat)/(r_sta+r_sat-r_sta2sat))
            Rela=Rela - dbending
        end if
    end if
 
    ! Calculate the satellite clock correction
    i=1
    do while (i<=2)
        j=minloc(dabs(TT),dim=1)
       Clk(i)=SP3Data%Eph(PRN)%Clk(j)
        TTT(i)=TT(j)
          if (dabs(Clk(i)-999999.999999)>1.d0 )then
              i=i+1
              TT(j)=999999.d0
          else
              Sat_Clk=9999.d0
              return
          end if
    end do
    Sat_Clk=Clk(1)+(-t1-TTT(1))/(TTT(2)-TTT(1))*(Clk(2)-Clk(1))
    Sat_Clk=Sat_Clk*1.d-6
end subroutine