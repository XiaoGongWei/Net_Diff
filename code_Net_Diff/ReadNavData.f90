! =============  Read GPS Navigation Data  ==========
!
! PURPOSE:
!    Read GPS Navigation data from navidation ephemeris file.
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!  ===============  End of Header  =============

subroutine ReadNavData
use MOD_FileID
use MOD_NavData
use MOD_NavHead
use MOD_Var
use MOD_constant
implicit none
    integer :: Version
    character(80) line
    integer ::PRN, year, mon, day, hour, min, GPSweek, GPSweek0, Flag(32), Flag_J(7), Flag_I(7)
    real(8) :: sec, GPSsec, GPSsec0, dGPST
    
!    NavData(:,1)=NavData(:,2)
!    GPSweek0=3000
!    GPSsec0=0.d0
    Flag=1
    Flag_J=1
    Flag_I=1
    if (NavHead%Version==2) then
        do while(.true.)
            read(NavID,fmt="(A)",end=200) line  ! Line1
                read(line,"(I2,5I3,F5.1)") PRN, year, mon, day, hour, min, sec
                call UTC2GPST(year,mon, day, hour, min, sec, GPSweek, GPSsec)
                
                Flag(PRN)=Flag(PRN)+1
                if ( dabs((GPSweek-NavData(PRN)%Nav(Flag(PRN)-1)%GPSweek)*604800.0d0+GPSsec-NavData(PRN)%Nav(Flag(PRN)-1)%GPSsec)<0.1d0 .and. &
                        NavData(PRN)%Nav(Flag(PRN)-1)%Health/=0.d0 ) then
                    Flag(PRN)=Flag(PRN)-1 ! cover last unhealthy ephemeris of the same epoch
                end if 
                NavData(PRN)%Nav(Flag(PRN))%GPSweek=GPSweek
                NavData(PRN)%Nav(Flag(PRN))%GPSsec=GPSsec
                read(line,"(22X,3D19.12)") NavData(PRN)%Nav(Flag(PRN))%a0, NavData(PRN)%Nav(Flag(PRN))%a1, NavData(PRN)%Nav(Flag(PRN))%a2
            read(NavID,fmt="(A)",end=200) line ! Line2
                read(line,"(3X,19X,3D19.12)") NavData(PRN)%Nav(Flag(PRN))%Crs, NavData(PRN)%Nav(Flag(PRN))%delN, NavData(PRN)%Nav(Flag(PRN))%M0
            read(NavID,fmt="(A)",end=200) line ! Line3
                read(line,"(3X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%Cuc, NavData(PRN)%Nav(Flag(PRN))%e, NavData(PRN)%Nav(Flag(PRN))%Cus, NavData(PRN)%Nav(Flag(PRN))%sqrtA
            read(NavID,fmt="(A)",end=200) line ! Line4
                read(line,"(3X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%toe,NavData(PRN)%Nav(Flag(PRN))%Cic, NavData(PRN)%Nav(Flag(PRN))%Omega, NavData(PRN)%Nav(Flag(PRN))%Cis
            read(NavID,fmt="(A)",end=200) line ! Line5
                read(line,"(3X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%i0, NavData(PRN)%Nav(Flag(PRN))%Crc, NavData(PRN)%Nav(Flag(PRN))%w, NavData(PRN)%Nav(Flag(PRN))%Omegadot
            read(NavID,fmt="(A)",end=200) line ! Line6
                read(line,"(3X,D19.12,19X,D19.12)") NavData(PRN)%Nav(Flag(PRN))%idot, NavData(PRN)%Nav(Flag(PRN))%WeekNo
            read(NavID,fmt="(A)",end=200) line ! Line7
                read(line,"(3X,19X,2D19.12)") NavData(PRN)%Nav(Flag(PRN))%Health, NavData(PRN)%Nav(Flag(PRN))%TGD(1)
            read(NavID,fmt="(A)",end=200) line ! Line8
        end do
    elseif (NavHead%Version==3) then
        do while(.true.)
            read(NavID,fmt="(A)",end=200) line  ! Line1
            if ((line(1:2)=="GP") .or.  (line(1:2)=="GA")) cycle
            if ( (line(1:1)=="G") .and. (SystemUsed(1)) .and. (.not.(IsCNAV)) ) then
                read(line,"(1X,I2,I5,4I3,F3.0)")  PRN, year, mon, day, hour, min, sec
                if (PRN>GNum) cycle
                call UTC2GPST(year,mon, day, hour, min, sec, GPSweek, GPSsec)
                
                Flag(PRN)=Flag(PRN)+1
                if ( dabs((GPSweek-NavData(PRN)%Nav(Flag(PRN)-1)%GPSweek)*604800.0d0+GPSsec-NavData(PRN)%Nav(Flag(PRN)-1)%GPSsec)<0.1d0 .and. &
                        NavData(PRN)%Nav(Flag(PRN)-1)%Health/=0.d0 ) then
                    Flag(PRN)=Flag(PRN)-1 ! cover last unhealthy ephemeris of the same epoch
                end if 
                NavData(PRN)%Nav(Flag(PRN))%GPSweek=GPSweek
                NavData(PRN)%Nav(Flag(PRN))%GPSsec=GPSsec
                read(line,"(23X,3D19.12)") NavData(PRN)%Nav(Flag(PRN))%a0, NavData(PRN)%Nav(Flag(PRN))%a1, NavData(PRN)%Nav(Flag(PRN))%a2
                read(NavID,fmt="(A)",end=200) line ! Line2
                    read(line,"(4X,19X,3D19.12)") NavData(PRN)%Nav(Flag(PRN))%Crs, NavData(PRN)%Nav(Flag(PRN))%delN, NavData(PRN)%Nav(Flag(PRN))%M0
                read(NavID,fmt="(A)",end=200) line ! Line3
                    read(line,"(4X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%Cuc, NavData(PRN)%Nav(Flag(PRN))%e, NavData(PRN)%Nav(Flag(PRN))%Cus, NavData(PRN)%Nav(Flag(PRN))%sqrtA
                read(NavID,fmt="(A)",end=200) line ! Line4
                    read(line,"(4X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%toe,NavData(PRN)%Nav(Flag(PRN))%Cic, NavData(PRN)%Nav(Flag(PRN))%Omega, NavData(PRN)%Nav(Flag(PRN))%Cis
                read(NavID,fmt="(A)",end=200) line ! Line5
                    read(line,"(4X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%i0, NavData(PRN)%Nav(Flag(PRN))%Crc, NavData(PRN)%Nav(Flag(PRN))%w, NavData(PRN)%Nav(Flag(PRN))%Omegadot
                read(NavID,fmt="(A)",end=200) line ! Line6
                    read(line,"(4X,D19.12,19X,D19.12)") NavData(PRN)%Nav(Flag(PRN))%idot, NavData(PRN)%Nav(Flag(PRN))%WeekNo
                read(NavID,fmt="(A)",end=200) line ! Line7
                    read(line,"(4X,19X,2D19.12)") NavData(PRN)%Nav(Flag(PRN))%Health, NavData(PRN)%Nav(Flag(PRN))%TGD(1)
                read(NavID,fmt="(A)",end=200) line ! Line8
            elseif ( (line(1:1)=="J") .and. (SystemUsed(5)) ) then  ! QZSS
                read(line,"(1X,I2,I5,4I3,F3.0)")  PRN, year, mon, day, hour, min, sec
                if (PRN>JNum) cycle
                call UTC2GPST(year,mon, day, hour, min, sec, GPSweek, GPSsec)
                
                Flag_J(PRN)=Flag_J(PRN)+1
                if ( dabs((GPSweek-NavData_J(PRN)%Nav(Flag_J(PRN)-1)%GPSweek)*604800.0d0+GPSsec-NavData_J(PRN)%Nav(Flag_J(PRN)-1)%GPSsec)<0.1d0 .and. &
                        NavData_J(PRN)%Nav(Flag_J(PRN)-1)%Health/=0.d0 ) then
                    Flag_J(PRN)=Flag_J(PRN)-1 ! cover last unhealthy ephemeris of the same epoch
                end if 
                NavData_J(PRN)%Nav(Flag_J(PRN))%GPSweek=GPSweek
                NavData_J(PRN)%Nav(Flag_J(PRN))%GPSsec=GPSsec
                read(line,"(23X,3D19.12)") NavData_J(PRN)%Nav(Flag_J(PRN))%a0, NavData_J(PRN)%Nav(Flag_J(PRN))%a1, NavData_J(PRN)%Nav(Flag_J(PRN))%a2
                read(NavID,fmt="(A)",end=200) line ! Line2
                    read(line,"(4X,19X,3D19.12)") NavData_J(PRN)%Nav(Flag_J(PRN))%Crs, NavData_J(PRN)%Nav(Flag_J(PRN))%delN, NavData_J(PRN)%Nav(Flag_J(PRN))%M0
                read(NavID,fmt="(A)",end=200) line ! Line3
                    read(line,"(4X,4D19.12)") NavData_J(PRN)%Nav(Flag_J(PRN))%Cuc, NavData_J(PRN)%Nav(Flag_J(PRN))%e, NavData_J(PRN)%Nav(Flag_J(PRN))%Cus, NavData_J(PRN)%Nav(Flag_J(PRN))%sqrtA
                read(NavID,fmt="(A)",end=200) line ! Line4
                    read(line,"(4X,4D19.12)") NavData_J(PRN)%Nav(Flag_J(PRN))%toe,NavData_J(PRN)%Nav(Flag_J(PRN))%Cic, NavData_J(PRN)%Nav(Flag_J(PRN))%Omega, NavData_J(PRN)%Nav(Flag_J(PRN))%Cis
                read(NavID,fmt="(A)",end=200) line ! Line5
                    read(line,"(4X,4D19.12)") NavData_J(PRN)%Nav(Flag_J(PRN))%i0, NavData_J(PRN)%Nav(Flag_J(PRN))%Crc, NavData_J(PRN)%Nav(Flag_J(PRN))%w, NavData_J(PRN)%Nav(Flag_J(PRN))%Omegadot
                read(NavID,fmt="(A)",end=200) line ! Line6
                    read(line,"(4X,D19.12,19X,D19.12)") NavData_J(PRN)%Nav(Flag_J(PRN))%idot, NavData_J(PRN)%Nav(Flag_J(PRN))%WeekNo
                read(NavID,fmt="(A)",end=200) line ! Line7
                    read(line,"(4X,19X,2D19.12)") NavData_J(PRN)%Nav(Flag_J(PRN))%Health, NavData_J(PRN)%Nav(Flag_J(PRN))%TGD(1)
                read(NavID,fmt="(A)",end=200) line ! Line8
                if ( NavData_J(PRN)%Nav(Flag_J(PRN))%WeekNo==0.d0 .or. NavData_J(PRN)%Nav(Flag_J(PRN))%e==0.d0 .or. NavData_J(PRN)%Nav(Flag_J(PRN))%sqrtA==0.d0 ) then
                    Flag_J(PRN)=Flag_J(PRN)-1
                end if 
            elseif ( (line(1:1)=="I") .and. (SystemUsed(6)) ) then  ! IRNSS
                read(line,"(1X,I2,I5,4I3,F3.0)")  PRN, year, mon, day, hour, min, sec
                if (PRN>INum) cycle
                call UTC2GPST(year,mon, day, hour, min, sec, GPSweek, GPSsec)
                
                Flag_I(PRN)=Flag_I(PRN)+1
                if ( dabs((GPSweek-NavData_I(PRN)%Nav(Flag_I(PRN)-1)%GPSweek)*604800.0d0+GPSsec-NavData_I(PRN)%Nav(Flag_I(PRN)-1)%GPSsec)<0.1d0 .and. &
                        NavData_I(PRN)%Nav(Flag_I(PRN)-1)%Health/=0.d0 ) then
                    Flag_I(PRN)=Flag_I(PRN)-1 ! cover last unhealthy ephemeris of the same epoch
                end if 
                NavData_I(PRN)%Nav(Flag_I(PRN))%GPSweek=GPSweek
                NavData_I(PRN)%Nav(Flag_I(PRN))%GPSsec=GPSsec
                read(line,"(23X,3D19.12)") NavData_I(PRN)%Nav(Flag_I(PRN))%a0, NavData_I(PRN)%Nav(Flag_I(PRN))%a1, NavData_I(PRN)%Nav(Flag_I(PRN))%a2
                read(NavID,fmt="(A)",end=200) line ! Line2
                    read(line,"(4X,19X,3D19.12)") NavData_I(PRN)%Nav(Flag_I(PRN))%Crs, NavData_I(PRN)%Nav(Flag_I(PRN))%delN, NavData_I(PRN)%Nav(Flag_I(PRN))%M0
                read(NavID,fmt="(A)",end=200) line ! Line3
                    read(line,"(4X,4D19.12)") NavData_I(PRN)%Nav(Flag_I(PRN))%Cuc, NavData_I(PRN)%Nav(Flag_I(PRN))%e, NavData_I(PRN)%Nav(Flag_I(PRN))%Cus, NavData_I(PRN)%Nav(Flag_I(PRN))%sqrtA
                read(NavID,fmt="(A)",end=200) line ! Line4
                    read(line,"(4X,4D19.12)") NavData_I(PRN)%Nav(Flag_I(PRN))%toe,NavData_I(PRN)%Nav(Flag_I(PRN))%Cic, NavData_I(PRN)%Nav(Flag_I(PRN))%Omega, NavData_I(PRN)%Nav(Flag_I(PRN))%Cis
                read(NavID,fmt="(A)",end=200) line ! Line5
                    read(line,"(4X,4D19.12)") NavData_I(PRN)%Nav(Flag_I(PRN))%i0, NavData_I(PRN)%Nav(Flag_I(PRN))%Crc, NavData_I(PRN)%Nav(Flag_I(PRN))%w, NavData_I(PRN)%Nav(Flag_I(PRN))%Omegadot
                read(NavID,fmt="(A)",end=200) line ! Line6
                    read(line,"(4X,D19.12,19X,D19.12)") NavData_I(PRN)%Nav(Flag_I(PRN))%idot, NavData_I(PRN)%Nav(Flag_I(PRN))%WeekNo
                read(NavID,fmt="(A)",end=200) line ! Line7
                    read(line,"(4X,19X,2D19.12)") NavData_I(PRN)%Nav(Flag_I(PRN))%Health, NavData_I(PRN)%Nav(Flag_I(PRN))%TGD(1)
                read(NavID,fmt="(A)",end=200) line ! Line8
                if ( NavData_I(PRN)%Nav(Flag_I(PRN))%WeekNo==0.d0 .or. NavData_I(PRN)%Nav(Flag_I(PRN))%e==0.d0 .or. NavData_I(PRN)%Nav(Flag_I(PRN))%sqrtA==0.d0 ) then
                    Flag_I(PRN)=Flag_I(PRN)-1
                end if 
            end if
        end do
    elseif (NavHead%Version==4) then
        do while(.true.)
            read(NavID,fmt="(A)",end=200) line  ! Line1
            if ( (line(1:1)=='G') .and. (SystemUsed(1)) .and. (index(line,'CNAV')/=0) ) then
                read(line,"(1X,I2)")  PRN
                if (PRN>GNum) cycle
                read(NavID,fmt="(A)",end=200) line  ! Line1
                read(line,"(3X,I5,4I3,F3.0)")  year, mon, day, hour, min, sec
                call UTC2GPST(year,mon, day, hour, min, sec, GPSweek, GPSsec)
                
                Flag(PRN)=Flag(PRN)+1
                if ( dabs((GPSweek-NavData(PRN)%Nav(Flag(PRN)-1)%GPSweek)*604800.0d0+GPSsec-NavData(PRN)%Nav(Flag(PRN)-1)%GPSsec)<0.1d0 .and. &
                        NavData(PRN)%Nav(Flag(PRN)-1)%Health/=0.d0 ) then
                    Flag(PRN)=Flag(PRN)-1 ! cover last unhealthy ephemeris of the same epoch
                end if 
                NavData(PRN)%Nav(Flag(PRN))%GPSweek=GPSweek
                NavData(PRN)%Nav(Flag(PRN))%GPSsec=GPSsec
                read(line,"(23X,3D19.12)") NavData(PRN)%Nav(Flag(PRN))%a0, NavData(PRN)%Nav(Flag(PRN))%a1, NavData(PRN)%Nav(Flag(PRN))%a2
                read(NavID,fmt="(A)",end=200) line ! Line2
                    read(line,"(4X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%aDot, NavData(PRN)%Nav(Flag(PRN))%Crs, NavData(PRN)%Nav(Flag(PRN))%delN, NavData(PRN)%Nav(Flag(PRN))%M0
                read(NavID,fmt="(A)",end=200) line ! Line3
                    read(line,"(4X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%Cuc, NavData(PRN)%Nav(Flag(PRN))%e, NavData(PRN)%Nav(Flag(PRN))%Cus, NavData(PRN)%Nav(Flag(PRN))%sqrtA
                read(NavID,fmt="(A)",end=200) line ! Line4
                    read(line,"(4X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%toe,NavData(PRN)%Nav(Flag(PRN))%Cic, NavData(PRN)%Nav(Flag(PRN))%Omega, NavData(PRN)%Nav(Flag(PRN))%Cis
                    NavData(PRN)%Nav(Flag(PRN))%toe=GPSsec  ! different from LNAV
                read(NavID,fmt="(A)",end=200) line ! Line5
                    read(line,"(4X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%i0, NavData(PRN)%Nav(Flag(PRN))%Crc, NavData(PRN)%Nav(Flag(PRN))%w, NavData(PRN)%Nav(Flag(PRN))%Omegadot
                read(NavID,fmt="(A)",end=200) line ! Line6
                    read(line,"(4X,2D19.12)") NavData(PRN)%Nav(Flag(PRN))%idot, NavData(PRN)%Nav(Flag(PRN))%ndot
                read(NavID,fmt="(A)",end=200) line ! Line7
                    read(line,"(4X,19X,2D19.12)") NavData(PRN)%Nav(Flag(PRN))%Health, NavData(PRN)%Nav(Flag(PRN))%TGD(1)
                read(NavID,fmt="(A)",end=200) line ! Line8
                    read(line,"(4X,4D19.12)") NavData(PRN)%Nav(Flag(PRN))%ISCL1CA, NavData(PRN)%Nav(Flag(PRN))%ISCL2C, NavData(PRN)%Nav(Flag(PRN))%ISCL5I5, NavData(PRN)%Nav(Flag(PRN))%ISCL5Q5
                read(NavID,fmt="(A)",end=200) line ! Line9
            end if
        end do
    end if
    200 close(NavID)
    return
end subroutine