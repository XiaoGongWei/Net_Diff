! ==========Read GLONASS Navigation Data========
!
! PURPOSE:
!    Read GPS Navigation data from navidation ephemeris file.
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!  ===============  End of Header  =============

subroutine ReadNavData_R
use MOD_FileID
use MOD_NavData_R
use MOD_NavHead
use MOD_constant
implicit none
    character(80) line
    integer ::PRN, year, mon, day, hour, min, GPSweek, GPSweek0, Flag(27)
    real(8) :: sec, GPSsec, GPSsec0, dGPST
    
!    NavData(:,1)=NavData(:,2)
!    GPSweek0=3000
!    GPSsec0=0.d0
    Flag=1
    if (NavHead_R%Version==2) then
        do while(.true.)
            read(NavID_R,fmt="(A)",end=200) line  ! Line1
                read(line,"(I2,5I3,F5.1)") PRN, year, mon, day, hour, min, sec
                call UTC2GPST(year,mon, day, hour, min, sec, GPSweek, GPSsec)
                
                dGPST=(GPSweek-NavData_R(PRN)%Nav(Flag(PRN))%GPSweek)*604800.d0+GPSsec-NavData_R(PRN)%Nav(Flag(PRN))%GPSsec   ! Time difference between two PRN
    !            if (dGPST>1000.d0) then
                    Flag(PRN)=Flag(PRN)+1
    !            end if
                NavData_R(PRN)%Nav(Flag(PRN))%year=year
                NavData_R(PRN)%Nav(Flag(PRN))%mon=mon
                NavData_R(PRN)%Nav(Flag(PRN))%day=day
                NavData_R(PRN)%Nav(Flag(PRN))%hour=hour
                NavData_R(PRN)%Nav(Flag(PRN))%min=min
                NavData_R(PRN)%Nav(Flag(PRN))%sec=sec
                NavData_R(PRN)%Nav(Flag(PRN))%GPSweek=GPSweek
                NavData_R(PRN)%Nav(Flag(PRN))%GPSsec=GPSsec
                read(line,"(22X,2D19.12)") NavData_R(PRN)%Nav(Flag(PRN))%a0, NavData_R(PRN)%Nav(Flag(PRN))%a1
            read(NavID_R,fmt="(A)",end=200) line ! Line2
                read(line,"(3X,4D19.12)") NavData_R(PRN)%Nav(Flag(PRN))%X, NavData_R(PRN)%Nav(Flag(PRN))%Vx, NavData_R(PRN)%Nav(Flag(PRN))%Ax, NavData_R(PRN)%Nav(Flag(PRN))%Health
            read(NavID_R,fmt="(A)",end=200) line ! Line3
                read(line,"(3X,3D19.12)") NavData_R(PRN)%Nav(Flag(PRN))%Y, NavData_R(PRN)%Nav(Flag(PRN))%Vy, NavData_R(PRN)%Nav(Flag(PRN))%Ay
            read(NavID_R,fmt="(A)",end=200) line ! Line4
                read(line,"(3X,3D19.12)") NavData_R(PRN)%Nav(Flag(PRN))%Z,NavData_R(PRN)%Nav(Flag(PRN))%Vz, NavData_R(PRN)%Nav(Flag(PRN))%Az
        end do
    elseif (NavHead_R%Version==3) then
        do while(.true.)
            read(NavID_R,fmt="(A)",end=200) line  ! Line1
            if (line(1:1)=="R") then
                read(line,"(1X,I2,I5,4I3,F3.0)")  PRN, year, mon, day, hour, min, sec
                if (PRN>RNum) cycle
                call UTC2GPST(year,mon, day, hour, min, sec, GPSweek, GPSsec)
                
                dGPST=(GPSweek-NavData_R(PRN)%Nav(Flag(PRN))%GPSweek)*604800.d0+GPSsec-NavData_R(PRN)%Nav(Flag(PRN))%GPSsec   ! Time difference between two PRN
    !            if (dGPST>1000.d0) then
                    Flag(PRN)=Flag(PRN)+1
    !            end if
                NavData_R(PRN)%Nav(Flag(PRN))%year=year
                NavData_R(PRN)%Nav(Flag(PRN))%mon=mon
                NavData_R(PRN)%Nav(Flag(PRN))%day=day
                NavData_R(PRN)%Nav(Flag(PRN))%hour=hour
                NavData_R(PRN)%Nav(Flag(PRN))%min=min
                NavData_R(PRN)%Nav(Flag(PRN))%sec=sec
                NavData_R(PRN)%Nav(Flag(PRN))%GPSweek=GPSweek
                NavData_R(PRN)%Nav(Flag(PRN))%GPSsec=GPSsec
                read(line,"(23X,2D19.12)") NavData_R(PRN)%Nav(Flag(PRN))%a0, NavData_R(PRN)%Nav(Flag(PRN))%a1
                read(NavID_R,fmt="(A)",end=200) line ! Line2
                read(line,"(4X,4D19.12)") NavData_R(PRN)%Nav(Flag(PRN))%X, NavData_R(PRN)%Nav(Flag(PRN))%Vx, NavData_R(PRN)%Nav(Flag(PRN))%Ax, NavData_R(PRN)%Nav(Flag(PRN))%Health
                read(NavID_R,fmt="(A)",end=200) line ! Line3
                read(line,"(4X,3D19.12)") NavData_R(PRN)%Nav(Flag(PRN))%Y, NavData_R(PRN)%Nav(Flag(PRN))%Vy, NavData_R(PRN)%Nav(Flag(PRN))%Ay
                read(NavID_R,fmt="(A)",end=200) line ! Line4
                read(line,"(4X,3D19.12)") NavData_R(PRN)%Nav(Flag(PRN))%Z,NavData_R(PRN)%Nav(Flag(PRN))%Vz, NavData_R(PRN)%Nav(Flag(PRN))%Az
            end if

        end do
    end if
    200 close(NavID_R)
    return
end subroutine