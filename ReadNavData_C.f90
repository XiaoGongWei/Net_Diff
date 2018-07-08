! ==========  Read BEIDOU Navigation Data  ========
!
! PURPOSE:
!    Read BEIDOU Navigation data from navidation ephemeris file.
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!  ===============  End of Header  =============

subroutine ReadNavData_C
use MOD_FileID
use MOD_NavData
use MOD_NavHead
use MOD_constant
implicit none
    character(80) :: line
    integer ::PRN, year, mon, day, hour, min, GPSweek, GPSweek0, Flag(35)
    real(8) :: sec, GPSsec, GPSsec0, dGPST
    
!    NavData(:,1)=NavData(:,2)
!    GPSweek0=3000
!    GPSsec0=0.d0
    Flag=1
    if (NavHead_C%Version==2) then
        do while(.true.)
            read(NavID_C,fmt="(A)",end=200) line  ! Line1
                read(line,"(I2,5I3,F5.1)") PRN, year, mon, day, hour, min, sec
                call UTC2GPST(year,mon, day, hour, min, sec, GPSweek, GPSsec)

                Flag(PRN)=Flag(PRN)+1
                if ( dabs((GPSweek-NavData_C(PRN)%Nav(Flag(PRN)-1)%GPSweek)*604800.0d0+GPSsec-NavData_C(PRN)%Nav(Flag(PRN)-1)%GPSsec)<0.1d0 .and. &
                        NavData_C(PRN)%Nav(Flag(PRN)-1)%Health/=0.d0 ) then
                    Flag(PRN)=Flag(PRN)-1 ! cover last unhealthy ephemeris of the same epoch
                end if 
                NavData_C(PRN)%Nav(Flag(PRN))%GPSweek=GPSweek
                NavData_C(PRN)%Nav(Flag(PRN))%GPSsec=GPSsec
                read(line,"(22X,3D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%a0, NavData_C(PRN)%Nav(Flag(PRN))%a1, NavData_C(PRN)%Nav(Flag(PRN))%a2
                
            read(NavID_C,fmt="(A)",end=200) line ! Line2
                read(line,"(3X,4D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%IODE, NavData_C(PRN)%Nav(Flag(PRN))%Crs, NavData_C(PRN)%Nav(Flag(PRN))%delN, NavData_C(PRN)%Nav(Flag(PRN))%M0
            read(NavID_C,fmt="(A)",end=200) line ! Line3
                read(line,"(3X,4D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%Cuc, NavData_C(PRN)%Nav(Flag(PRN))%e, NavData_C(PRN)%Nav(Flag(PRN))%Cus, NavData_C(PRN)%Nav(Flag(PRN))%sqrtA
            read(NavID_C,fmt="(A)",end=200) line ! Line4
                read(line,"(3X,4D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%toe,NavData_C(PRN)%Nav(Flag(PRN))%Cic, NavData_C(PRN)%Nav(Flag(PRN))%Omega, NavData_C(PRN)%Nav(Flag(PRN))%Cis
            read(NavID_C,fmt="(A)",end=200) line ! Line5
                read(line,"(3X,4D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%i0, NavData_C(PRN)%Nav(Flag(PRN))%Crc, NavData_C(PRN)%Nav(Flag(PRN))%w, NavData_C(PRN)%Nav(Flag(PRN))%Omegadot
            read(NavID_C,fmt="(A)",end=200) line ! Line6
                read(line,"(3X,D19.12,19X,D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%idot, NavData_C(PRN)%Nav(Flag(PRN))%WeekNo
            read(NavID_C,fmt="(A)",end=200) line ! Line7
                read(line,"(3X,19X,3D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%Health, NavData_C(PRN)%Nav(Flag(PRN))%TGD(1), NavData_C(PRN)%Nav(Flag(PRN))%TGD(2)
            read(NavID_C,fmt="(A)",end=200) line ! Line8
                read(line,"(3X,19X,D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%IODC
            if ( NavData_C(PRN)%Nav(Flag(PRN))%WeekNo==0.d0 .or. NavData_C(PRN)%Nav(Flag(PRN))%e==0.d0 .or. NavData_C(PRN)%Nav(Flag(PRN))%sqrtA==0.d0 ) then
                Flag(PRN)=Flag(PRN)-1
            end if 
        end do
    elseif (NavHead_C%Version==3) then
        do while(.true.)
            read(NavID_C,fmt="(A)",end=200) line  ! Line1
            if (line(1:1)=="C") then
                read(line,"(1X,I2,I5,4I3,F3.0)")  PRN, year, mon, day, hour, min, sec
                if (PRN>CNum) cycle
                !sec=sec+14.d0   ! 广播星历时间为BDT，观测时间为GPST
                call UTC2GPST(year,mon, day, hour, min, sec, GPSweek, GPSsec)

                Flag(PRN)=Flag(PRN)+1
                if ( dabs((GPSweek-NavData_C(PRN)%Nav(Flag(PRN)-1)%GPSweek)*604800.0d0+GPSsec-NavData_C(PRN)%Nav(Flag(PRN)-1)%GPSsec)<0.1d0 .and. &
                    NavData_C(PRN)%Nav(Flag(PRN)-1)%Health/=0.d0 ) then
                    Flag(PRN)=Flag(PRN)-1 ! cover last unhealthy ephemeris of the same epoch
                end if 
                NavData_C(PRN)%Nav(Flag(PRN))%GPSweek=GPSweek
                NavData_C(PRN)%Nav(Flag(PRN))%GPSsec=GPSsec
                read(line,"(23X,3D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%a0, NavData_C(PRN)%Nav(Flag(PRN))%a1, NavData_C(PRN)%Nav(Flag(PRN))%a2
                read(NavID_C,fmt="(A)",end=200) line ! Line2
                    read(line,"(4X,4D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%IODE, NavData_C(PRN)%Nav(Flag(PRN))%Crs, NavData_C(PRN)%Nav(Flag(PRN))%delN, NavData_C(PRN)%Nav(Flag(PRN))%M0
                read(NavID_C,fmt="(A)",end=200) line ! Line3
                    read(line,"(4X,4D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%Cuc, NavData_C(PRN)%Nav(Flag(PRN))%e, NavData_C(PRN)%Nav(Flag(PRN))%Cus, NavData_C(PRN)%Nav(Flag(PRN))%sqrtA
                read(NavID_C,fmt="(A)",end=200) line ! Line4
                    read(line,"(4X,4D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%toe,NavData_C(PRN)%Nav(Flag(PRN))%Cic, NavData_C(PRN)%Nav(Flag(PRN))%Omega, NavData_C(PRN)%Nav(Flag(PRN))%Cis
                read(NavID_C,fmt="(A)",end=200) line ! Line5
                    read(line,"(4X,4D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%i0, NavData_C(PRN)%Nav(Flag(PRN))%Crc, NavData_C(PRN)%Nav(Flag(PRN))%w, NavData_C(PRN)%Nav(Flag(PRN))%Omegadot
                read(NavID_C,fmt="(A)",end=200) line ! Line6
                    read(line,"(4X,D19.12,19X,D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%idot, NavData_C(PRN)%Nav(Flag(PRN))%WeekNo
                read(NavID_C,fmt="(A)",end=200) line ! Line7
                    read(line,"(4X,19X,3D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%Health, NavData_C(PRN)%Nav(Flag(PRN))%TGD(1),  NavData_C(PRN)%Nav(Flag(PRN))%TGD(2)
                read(NavID_C,fmt="(A)",end=200) line ! Line8
                    read(line,"(4X,19X,D19.12)") NavData_C(PRN)%Nav(Flag(PRN))%IODC
                if ( NavData_C(PRN)%Nav(Flag(PRN))%WeekNo==0.d0 .or.  NavData_C(PRN)%Nav(Flag(PRN))%e==0.d0 .or. NavData_C(PRN)%Nav(Flag(PRN))%sqrtA==0.d0) then
                    Flag(PRN)=Flag(PRN)-1
                end if 
            end if    ! if (line(1:1)=="C") then
        end do
    end if   !  if (NavHead_C%Version==2) then
    200 close(NavID_C)
    return
end subroutine