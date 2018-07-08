! ======================Cal_Sat_Pos_e====================
!
! PURPOSE:
!           Calculate the GALILEO Satellite Coordinate, satellite clock error, 
!     relativity error, the range between the satellite and the receiver, 
!     using GPS broadcast ephemeris.
!
! REFERENCE:
!    http://www.navipedia.net/index.php/Emission_Time_Computation
!    http://www.navipedia.net/index.php/GPS_and_Galileo_Satellite_Coordinates_Computation
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
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================End of Header==============

subroutine Cal_Sat_Pos_e(GPSweek, GPSsec, PRN, Rec_Coor, t1, Sat_Coor2, Sat_Vel, R, Rela)
use MOD_constant
use MOD_NavData
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
    
    type(type_NavPRN) :: tempNav
    real(8) :: t, tkr, tk
    real(8) :: e, A, n0,n
    real(8) :: Mk, Ek0, Ek, vk, PHI_k
    real(8) :: duk, dik, drk
    real(8) :: uk, radius, ik
    real(8) :: x, y, Omega
    real(8) :: F, Rotation(3,3), t2
    integer :: i
    real(8) :: EkDot, PHI_Dot, ukDot, rkDot, ikDot, xDot, yDot, zDot, OmegaDot

    tkr=(GPSweek-Navdata_E(PRN)%Nav(1)%GPSweek)*604800.d0+GPSsec-Navdata_E(PRN)%Nav(1)%GPSsec
    do i=2,Ubound(Navdata_E(PRN)%Nav,dim=1)
        t=(GPSweek-Navdata_E(PRN)%Nav(i)%GPSweek)*604800.d0+GPSsec-Navdata_E(PRN)%Nav(i)%GPSsec
!        if (dabs(t)<=dabs(tkr)-0.1d0) then
!            tkr=t    ! 取距离最近时刻的星历
!            tempNav=Navdata_E(PRN)%Nav(i)
!        else if ( (dabs(t) - dabs(tkr))>700.d0 ) then
!            exit
!        end if
        if ((t>=-0.1d0)  .and. (abs(t)-abs(tkr)<20.d0) .and. (NavData_E(PRN)%Nav(i)%Health==0.d0)) then
            tkr=t   ! 取最新的星历,t扣除了接收机钟差，!   .and. (mod(NavData_E(PRN)%Nav(i)%Code,2.d0)==0.d0) Code=258 means F/NAV, E1/E5a; Code =513 or 517 mens I/Nav, E1/E5b
            tempNav=Navdata_E(PRN)%Nav(i)
            if ((dabs(tkr)<3620.d0)  .and. (tempNav%Health==0.d0)  .and. (tempNav%a0/=0.d0) ) then
                exit
            end if
!        elseif (t<-0.2d0) then
!            exit
        end if
    end do
    if ( (dabs(tkr)>3620.d0) ) then
        Sat_Coor2=9999.d0
        Sat_Clk=9999.d0
        return
    end if
    omg=7.2921151467d-5        !  GTRF: Earth rotation angular velocity (rad/sec)
    GM=3.986004418d14   !  GTRF: universal gravitational param (m^3/s^2) ! See Galileo ICD
    e=tempNav%e
    A=tempNav%sqrtA**2
    n0=dsqrt(GM/A/A/A)
    n=n0+tempNav%delN  ! Corrected mean motion (rad/s)
    F= -4.442807309d-10   ! See Galileo ICD
    do while (.true.)
        tk=tkr-t1
        Mk=tempNav%M0+n*tk  ! Mean anomaly
        Ek0=Mk
        do while(.true.)
            Ek=Mk+e*dsin(Ek0)
            if (dabs(Ek-Ek0)<1.d-12) exit       ! Kepler's Equation for Eccentric Anomoly Ek
            Ek0=Ek
        end do
        EkDot=n/(1.d0-e*dcos(Ek))  ! For Satellite Velocity
        vk=2*datan(dtan(Ek/2.d0)*dsqrt((1.d0+e)/(1.d0-e)) )   ! True anomaly (rad)
        ! vk = atan2( sqrt(1.0-e*e)*sin(Ek), cos(Ek)-e);
        if (vk<0.d0) vk=vk+2*pi
        PHI_k=vk+tempNav%w   ! Argument of latitude
        PHI_Dot=sqrt(1.d0+e)/sqrt(1.d0-e)*(dcos(vk/2.d0))**2/(dcos(Ek/2.d0))**2*EkDot  ! For Satellite Velocity
        
        ! Second Harmonic Perturbations
        duk=tempNav%Cus*dsin(2*PHI_k)+tempNav%Cuc*dcos(2*PHI_k)   ! Argument of Lat correction
        drk=tempNav%Crs*dsin(2*PHI_k)+tempNav%Crc*dcos(2*PHI_k)     ! Radius correction
        dik=tempNav%Cis*dsin(2*PHI_k)+tempNav%Cic*dcos(2*PHI_k)      ! Inclination correction
        
        uk=PHI_k+duk                                                            ! Corr. arg of lat
        ukDot=(1.d0+2.d0*tempNav%Cus*dcos(2.d0*PHI_k)-2.d0*tempNav%Cuc*dsin(2*PHI_k))*PHI_Dot  ! For Satellite Velocity
        radius=A*(1-e*dcos(Ek))+drk                                             ! Corrected radius
        rkDot=EkDot*A*e*dsin(Ek)+2.d0*(tempNav%Crs*dcos(2.d0*PHI_k)-tempNav%Crc*dsin(2.d0*PHI_k))*PHI_Dot  ! For Satellite Velocity
        ik=tempNav%i0+dik+tempNav%idot*tk   ! Corrected inclination
        ikDot=2.d0*(tempNav%Cis*dcos(2.d0*PHI_k)-tempNav%Cic*dsin(2.d0*PHI_k))*PHI_Dot+tempNav%idot  ! For Satellite Velocity
        
        ! Positons in orbital plane
        x=radius*dcos(uk)
        y=radius*dsin(uk)
        xDot=rkDot*dcos(uk) - radius*dsin(uk)*ukDot  ! For Satellite Velocity
        yDot=rkDot*dsin(uk) + radius*dcos(uk)*ukDot  ! For Satellite Velocity
        
        Omega=tempNav%omega+(tempNav%Omegadot-omg)*tk-omg*tempNav%toe  ! Corrected longitude of ascending node
        OmegaDot=tempNav%Omegadot - omg  ! For Satellite Velocity
    
        ! ECEF coordinates
        Sat_Coor(1)=x*dcos(Omega)-y*dcos(ik)*dsin(Omega)
        Sat_Coor(2)=x*dsin(Omega)+y*dcos(ik)*dcos(Omega)
        Sat_Coor(3)=y*dsin(ik)
        Sat_Vel(1)=xDot*dcos(Omega)-yDot*dsin(Omega)*dcos(ik)+y*dsin(Omega)*dsin(ik)*ikDot &
                 - (x*dsin(Omega)+y*dcos(Omega)*dcos(ik))*OmegaDot  ! For Satellite Velocity
        Sat_Vel(2)=xDot*dsin(Omega)+yDot*dcos(Omega)*dcos(ik)-y*dcos(Omega)*dsin(ik)*ikDot &
                 + (x*dcos(Omega)-y*dsin(Omega)*dcos(ik))*OmegaDot  ! For Satellite Velocity
         Sat_Vel(3)=yDot*dsin(ik)+y*dcos(ik)*ikDot  ! For Satellite Velocity

        Rotation=reshape( (/ 0.d0, dsin(omg), 0, -dsin(omg),0.d0, 0.d0, 0.d0, 0.d0, 0.0d0 /), (/3,3/) ) 
        Sat_Vel=Sat_Vel - MATMUL(Sat_Coor,Rotation)
        
!        ! Satellite clock
!        Sat_Clk=tempNav%a0+tk*tempNav%a1+tk*tk*tempNav%a2
        
        ! Clock-error correction with theory of relativity
        Rela=c*F*e*dsqrt(A)*dsin(Ek)  ! Unit in meter  ! 见GPS ICD
        
        !  Earth rotation correction
        Rotation=reshape( (/ dcos(omg*t1), dsin(omg*t1), 0, -dsin(omg*t1),dcos(omg*t1), 0, 0, 0, 1.0d0 /), (/3,3/) ) 
        
        Sat_Coor2=MATMUL(Sat_Coor,Rotation)
        ! Diatance from the satellite to the reciver
        R=dsqrt(DOT_PRODUCT((Sat_Coor2-Rec_Coor),(Sat_Coor2-Rec_Coor)))
        ! Signal transfering time
        t2=R/c
        if (dabs(t2-t1)<1.d-10) exit
        t1=t2
    end do
end subroutine