! =================== Cal_Sat_Pos_g =================
!
! PURPOSE:
!           Calculate the GLONASS Satellite Coordinate, satellite clock error, 
!     relativity error, the range between the satellite and the receiver, 
!     using GLONASS broadcast ephemeris.
!
! REFERENCE:
!      http://www.navipedia.net/index.php/Emission_Time_Computation
!      http://www.navipedia.net/index.php/GLONASS_Satellite_Coordinates_Computation
!     http://gpsworld.com/directions-2014-new-horizons-of-glonass/
!    http://www.glonass-iac.ru/en/content/news/?ELEMENT_ID=721
!     GLONASS ICD 5.1(2008): https://www.unavco.org/help/glossary/docs/ICD_GLONASS_5.1_(2008)_en.pdf
!                   or http://www.spacecorp.ru/directions/glonass/control_document/
!
! NOTICE:
!       According to the GLONASS modernisation plan[3], the ephemeris information 
!   implementing the PZ-90.11 reference system was updated on all operational 
!   GLONASS satellites starting from 3:00 pm on December 31, 2013. From this 
!   time on, the satellites are broadcasting in the PZ-90.11[4]. This ECEF 
!   reference frame is an updated version of PZ-90, closest to the ITRF2000.
!       The transformation from PZ-90.11 to ITRF2008 contains only an origin shift 
!   vector, and no rotations nor scale factor
!        [X,Y,Z]ITRF2008=[X,Y,Z]PZ90.11+[0.003, 0.001, 0.02]m    Equ (2)
!
! INPUTS:
!    GPSweek           GPS week
!    GPSsec              GPS second
!    PRN                   satellite number, integer
!    Rec_Coor           coordinate of the station, unit in meter
!    t1                       signal transfer time using LC combined range, unit in second
!
! OUTPUTS:
!    Sat_Coor2          satellite coordinate, unit in meter
!    R                        range between the satellite and the receiver, unit in meter
!    Rela                   relativity correction, unit in meter
!    Sat_Clk              satellite clock, unit in second
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! =====================End of Header=================

subroutine Cal_Sat_Pos_g(GPSweek, GPSsec, PRN, Rec_Coor,t1, Sat_Coor2, Sat_Vel, R, Rela,Sat_Clk)
use MOD_constant
use MOD_NavData_R
use MOD_VAR
implicit none
    ! Intent in
    integer :: PRN
    integer :: GPSweek
    real(8) :: GPSsec
    real(8) :: Rec_Coor(3)
    real(8) ::  t1
    ! Intent out
    real(8) ::  Sat_Coor(3), R, Sat_Coor2(3)
    real(8) ::  Rela, Sat_Clk, Sat_Vel(3)
    
    real(8) :: tkr, t, tk, t2
    integer :: i, N
    type(type_NavPRN_R) :: tempNav
    real(8) :: X, Y, Z
    real(8) :: X0, Y0, Z0
    real(8) :: Vx, Vy, Vz
    real(8) :: Vx0, Vy0, Vz0
    real(8) :: Ax, Ay, Az
    real(8) :: Ax0, Ay0, Az0
    reaL(8) :: Step, Step1
    real(8) :: K11,K12,K13,K14,K15,K16
    real(8) :: K21,K22,K23,K24,K25,K26
    real(8) :: K31,K32,K33,K34,K35,K36
    real(8) :: K41,K42,K43,K44,K45,K46
    real(8) :: Rotation(3,3)
    
    tkr=(GPSweek-NavData_R(PRN)%Nav(1)%GPSweek)*604800.d0+GPSsec-NavData_R(PRN)%Nav(1)%GPSsec
    do i=2,Ubound(NavData_R(PRN)%Nav,dim=1)
        t=(GPSweek-NavData_R(PRN)%Nav(i)%GPSweek)*604800.d0+GPSsec-NavData_R(PRN)%Nav(i)%GPSsec
        if (dabs(t)<dabs(tkr) .and. (NavData_R(PRN)%Nav(i)%Health==0.d0)) then
            tkr=t   ! 取距离最近时刻的星历
            tempNav=NavData_R(PRN)%Nav(i)
        else if ( (dabs(t) - dabs(tkr))>700.d0 ) then
            exit
        end if
!        if ((t>=-0.1d0)  .and. (abs(t)-abs(tkr)<20.d0) .and. (NavData_R(PRN)%Nav(i)%Health==0.d0)) then 
!            tkr=t   ! 取最新的星历
!            tempNav=NavData_R(PRN)%Nav(i)
!        elseif (t<-0.2d0) then
!            exit
!        end if
    end do
    
   if ( (dabs(tkr)>1800.1d0) .or. (tempNav%X==0.d0) .or. (tempNav%Y==0.d0) .or. (tempNav%Z==0.d0) ) then
        Sat_Coor2=9999.d0
        Sat_Clk=9999.d0
        return
    end if
    
    do while (.true.)
        Step=30.d0
        tk=tkr-t1
        if(tk<0)then
            Step=-Step
        end if
        N=int(dabs(tk)/dabs(Step))
        X0=tempNav%X
        Y0=tempNav%Y
        Z0=tempNav%Z
        Vx0=tempNav%Vx
        Vy0=tempNav%Vy
        Vz0=tempNav%Vz
        Ax0=tempNav%Ax
        Ay0=tempNav%Ay
        Az0=tempNav%Az
        
        ! 4th order Runge-Kutta numerical integration 
        do i=0,N
            Step1=tk-i*Step
            if(dabs(Step1)<dabs(Step))then
                  Step=Step1
            endif
            ! K1: at X,Y
            Vx=Vx0
            Vy=Vy0
            Vz=Vz0
            X=X0
            Y=Y0
            Z=Z0
            call acceleration(X,Y,Z,Vx,Vy,Vz,Ax0,Ay0,Az0,K11,K12,K13)  ! acceleration at X
            K14=Vx
            K15=Vy
            K16=Vz     ! Velocity at X
            
            ! K2: at X+h/2, Y+K1*h/2
            Vx=Vx0+K11*Step/2.d0
            Vy=Vy0+K12*Step/2.d0
            Vz=Vz0+K13*Step/2.d0
!            X=X0+Vx*Step/2.d0   ! K14
!            Y=Y0+Vy*Step/2.d0    ! K15
!            Z=Z0+Vz*Step/2.d0    ! K16
            X=X0+K14*Step/2.d0   ! K14
            Y=Y0+K15*Step/2.d0    ! K15
            Z=Z0+K16*Step/2.d0    ! K16
            call acceleration(X,Y,Z,Vx,Vy,Vz,Ax0,Ay0,Az0,K21,K22,K23)  ! acceleration at X+Step/2
            K24=Vx
            K25=Vy
            K26=Vz     ! Velocity at X+Step/2
            
            ! K3: at X+h/2, Y+K2*h/2
            Vx=Vx0+K21*Step/2.d0
            Vy=Vy0+K22*Step/2.d0
            Vz=Vz0+K23*Step/2.d0
!            X=X0+Vx*Step/2.d0   !  K24
!            Y=Y0+Vy*Step/2.d0    ! K25
!            Z=Z0+Vz*Step/2.d0    ! K26
            X=X0+K24*Step/2.d0 
            Y=Y0+K25*Step/2.d0  
            Z=Z0+K26*Step/2.d0 
            call acceleration(X,Y,Z,Vx,Vy,Vz,Ax0,Ay0,Az0,K31,K32,K33)  ! acceleration at X+Step/2
            K34=Vx
            K35=Vy
            K36=Vz     ! Velocity at X+Step/2
            
            ! K4: at X+h, Y+K3*h
            Vx=Vx0+K31*Step
            Vy=Vy0+K32*Step
            Vz=Vz0+K33*Step
!            X=X0+Vx*Step     ! K34
!            Y=Y0+Vy*Step     ! K35
!            Z=Z0+Vz*Step    ! K36
            X=X0+K34*Step   
            Y=Y0+K35*Step 
            Z=Z0+K36*Step 
            call acceleration(X,Y,Z,Vx,Vy,Vz,Ax0,Ay0,Az0,K41,K42,K43)  ! acceleration at X+Step
            K44=Vx
            K45=Vy
            K46=Vz     ! Velocity at X+Step
            
            VX0=(K11+2.d0*K21+2.d0*K31+K41)*Step/6.d0+VX0
            VY0=(K12+2.d0*K22+2.d0*K32+K42)*Step/6.d0+VY0
            VZ0=(K13+2.d0*K23+2.d0*K33+K43)*Step/6.d0+VZ0
    
            X0=(K14+2.d0*K24+2.d0*K34+K44)*Step/6.d0+X0
            Y0=(K15+2.d0*K25+2.d0*K35+K45)*Step/6.d0+Y0
            Z0=(K16+2.d0*K26+2.d0*K36+K46)*Step/6.d0+Z0
        end do
        Sat_Coor=(/X, Y,Z/)*1000.d0
       
        !  Earth rotation correction
        Rotation=reshape( (/ dcos(omg*t1), dsin(omg*t1), 0, -dsin(omg*t1),dcos(omg*t1), 0, 0, 0, 1.0d0 /), (/3,3/) ) 
        Sat_Coor=MATMUL(Sat_Coor,Rotation)
        
        call PZ902WGS84(Sat_Coor, Sat_Coor2)   ! From PZ90 to WGS84
        R=dsqrt(DOT_PRODUCT((Sat_Coor2-Rec_Coor),(Sat_Coor2-Rec_Coor)))
        ! Signal transfering time
        t2=R/c
        if (dabs(t2-t1)<1.d-10) exit
        t1=t2
    end do
    
!    Rela=-(X*VX+Y*VY+Z*VZ)*2.d0/C*1000.d0*1000.d0
    Rela=0.d0
    Sat_Clk=tempNav%a0+tk*tempNav%a1

    
    ! TGD Correction
    ! GLONASS satelite clock is based on L1 C/A frequency
    ! Reference: Montenbruck, et al(2018). Multi-GNSS Signal-in-Space Range Error Assessment – Methodology and Results, Advance in Space Research
    if (index(ObsCombine,"PC")/=0)then
!        Sat_Clk=Sat_Clk -  f2**2/(f1+f2)/(f1-f2)*DCBBSX(4,PRN+GNum)   ! C1C-->C1C C2C
        Sat_Clk=Sat_Clk+ DCBBSX(3,PRN+GNum) -  f2**2/(f1+f2)/(f1-f2)*(DCBBSX(4,PRN+GNum)-DCBBSX(3,PRN+GNum))        !  C1C-->C1P--->C1P C2C
    elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) then
        Sat_Clk=Sat_Clk + DCBBSX(4,PRN+GNum)            !  C1C--->C2C
    elseif ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
        Sat_Clk=Sat_Clk + DCBBSX(3,PRN+GNum)            !  C1C-->C1P
    end if

    return
end subroutine


! ================ acceleration ======================
!
! PURPOSE:
!      Calculate the acceleration of the GLONASS satellite
!
! Reference:
!      http://www.navipedia.net/index.php/GLONASS_Satellite_Coordinates_Computation
! 
! INPUTS:
!      X, Y, Z               coordinate of the satellite
!      Vx, Vy, Vz          velocity of the satellite
!      Ax, Ay, Az         acceleration of the moon and sun 
!
! OUTPUTS:
!      Ax, Ay, Az          acceleration of the satellite
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================== End of Header =================

subroutine acceleration(X,Y,Z,Vx,Vy,Vz,Ax0,Ay0,Az0,Ax,Ay,Az)
implicit none
    ! Intent in:
    real(8) :: X, Y, Z
    real(8) :: Vx,Vy,Vz
    real(8) :: Ax0, Ay0, Az0
    ! Intent out:
    real(8) :: Ax, Ay, Az
    ! Local variables
    real(8), parameter :: u=398600.44d0
    real(8), parameter :: a=6378.136d0
    real(8), parameter :: C2= -1082.63d-6
    real(8), parameter :: omg=0.72921151d-4
    
    real(8) :: R,a1,a2
    
    R=dsqrt(X*X+Y*Y+Z*Z)
    a1=u/R/R/R
    a2=1.5d0*C2*u*a*a/R**5.d0
    Ax= -a1*X+a2*X*(1.d0-5*Z*Z/R/R)+omg*omg*X+2.d0*omg*Vy+Ax0
    Ay= -a1*Y+a2*Y*(1.d0-5*Z*Z/R/R)+omg*omg*Y-2.d0*omg*Vx+Ay0
    Az= -a1*Z+a2*Z*(3.d0-5*Z*Z/R/R)+Az0

    return
end subroutine


! =================== PZ902WGS84 ====================
!
! PURPOSE:
!         From PZ90.02 reference frame to WGS84.
! REFERENCES:
!     http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
!     http://gpsworld.com/directions-2014-new-horizons-of-glonass/
!    http://www.glonass-iac.ru/en/content/news/?ELEMENT_ID=721
!       GLONASS ICD 5.1(2008): https://www.unavco.org/help/glossary/docs/ICD_GLONASS_5.1_(2008)_en.pdf
!                   or http://www.spacecorp.ru/directions/glonass/control_document/
!     Vasty Engelsberg,Ivan Petrovski,Valery Babakov.
!     Glonass Business Prospects[J]. GPS WORLD,2008,(3):12-15.
!     任锴，杨力等.GPS/GLONASS组合单点定位精度分析[J].
!     海洋测绘,2012.30(1):11:13.
! 
! Notice:
!          According to the GLONASS modernisation plan, 
!    the ephemeris information implementing the PZ-90.02
!    reference system was updated on all operational GLONASS 
!    satellites from 12:00 to 17:00 UTC, September 20th., 2007. 
!    From this time on, the satellites are broadcasting in the 
!    PZ-90.02. This ECEF reference frame is an updated version 
!    of PZ-90, closest to the ITRF2000. 
!          The transformation from PZ-90.02 to ITRF2000 contains only an 
!   origin shift vector, but no rotations nor scale factor, as it is shown 
!   in equation (1) [Revnivykh, 2007]
!   [X,Y,Z]ITRF2000=[X,Y,Z]PZ90.02+[0.36, 0.08, 0.18]m    Equ (1)
!       [Revnivykh, 2007] Revnivykh, S., 2007. GLONASS Status and Progress. 
!   In: Minutes of the 47th CGSIC Meeting, Forth Worth, Texas.
!     According to the GLONASS modernisation plan[3], the ephemeris information 
!   implementing the PZ-90.11 reference system was updated on all operational 
!   GLONASS satellites starting from 3:00 pm on December 31, 2013. From this 
!   time on, the satellites are broadcasting in the PZ-90.11[4]. This ECEF 
!   reference frame is an updated version of PZ-90, closest to the ITRF2000.
!       The transformation from PZ-90.11 to ITRF2008 contains only an origin shift 
!   vector, and no rotations nor scale factor
!          [X,Y,Z]ITRF2008=[X,Y,Z]PZ90.11+[0.003, 0.001, 0.02]m    Equ (2)
!
! INPUT:
!     CoorPZ(3)              Coordinate in PZ90.11, unit in meter
!
! OUTPUT:
!     CoorWGS(3)          Coorinate in WGS84, unit in meter
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ======================End of Header======================

subroutine PZ902WGS84(CoorPZ,CoorWGS)
!use MOD_VAR
implicit none
! Intent in:
real(8) :: CoorPZ(3)
! Intent out:
real(8) :: CoorWGS(3)
! Local variables:
real(8) :: move(3)
    
!    if (int_day>=2014000) then
        move=(/0.003d0, 0.001d0, 0.02d0/)
!    elseif (int_day>=2007263) then
!        move=(/-0.36d0, 0.08d0, 0.18d0/)
!    else
!        write(*,*) "***GLONASS sysstem unknown****"
!    end if
    CoorWGS=CoorPZ+move
    return
end subroutine