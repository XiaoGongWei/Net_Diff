! ================  ReadObsData ==============
!
! PURPOSE:
!     Read observation data from **.**o file
!
! NOTICE:
!     1) The time system should be consistent with the
!     orbit and clock.
!     2) Both RINEX3.02 and 2.01 format can be handled here
!         For more information about RINEX file, refer to:
!            http://spot.colorado.edu/~kristine/rinex210.txt
!            ftp://igs.org/pub/data/format/rinex302.pdf
!
! INPUTS:
!      Obsweek         Obs time of GPS week
!      Obssec            Obs time of GPS second
!      sta                   order the multi-stations
!
! Output:
!      ObsData          Including C1, P1, P2, P3, L1, L2 data
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================  End of Header  ===============

subroutine ReadObsData(Obsweek, Obssec, ObsData,sta)
use MOD_ObsHead
use MOD_ObsData
use MOD_FileID
use MOD_VAR
implicit none
    type(type_ObsData) :: ObsData
    integer(1) :: Flag, lines
    integer :: Obsweek
    real(8) :: Obssec
    character*200 :: line
    character*800 :: line2
    integer :: year, mon, day, hour, min, GPSweek
    real(8) :: sec,int_sec, GPSsec
    integer :: i, j,k,sta, sys
    real(8) :: val
    character(1) :: System
    integer(1) :: PRN
    integer(1) :: LLI
    
    if (ObsHead(sta)%Version==2) then
        do while(.true.)
            ObsData%Flag=1   ! First set the flag false, in case that there is no data in this epoch
            read(unit=ObsID(sta),fmt="(A)",end=200) line
            if ( (index(line,"COMMENT") /=0) .or. (len_trim(line)==0) .or. (index(line,"SECONDS") /=0) ) then
                cycle
            end if
        
            read(line, "(28X,I1,I3)") Flag, lines
            if (Flag>1) then      ! Event flag
                do i=1,lines
                    read(unit=ObsID(sta), fmt="(A)",end=200) line
                    if (index(line(61:80),"ANTENNA: DELTA H/E/N") /=0 ) then  ! Event flag=3, eg. bjfs1660.10o
                        read(line,"(3F14.4)") ObsHead(sta)%Antenna(3), ObsHead(sta)%Antenna(2),ObsHead(sta)%Antenna(1)
                    end if
                end do
                cycle
            end if
        
            read(unit=line,fmt="(1X,I2,4(1X,I2),F11.7)") year, mon, day, hour, min, sec
            int_sec=real(nint(sec*10.d0))/10.d0   ! nint: nearest integer, for 1~10Hz data
            ObsData%Clk_Bias=sec-int_sec   ! station receiver clock correction
            ! int_sec=int_sec+14.d0 ! 精密星历时间为GPST，观测时间为BDT
        
            call UTC2GPST(year, mon, day, hour, min, int_sec, GPSweek, GPSsec)
        
            if ((GPSweek-Obsweek)*604800.d0+GPSsec-Obssec>0.01d0) then ! If no data in this epoch
            !if ((GPSweek>Obsweek) .or. (GPSsec>Obssec)) then ! If no data in this epoch
                backspace(ObsID(sta))
                exit
            end if
        
            if (year<80) then
                year=year+2000
            else if ((year>=80) .and. (year<100)) then
                year=year+1900
            end if
            ObsData%year=year
            ObsData%mon=mon
            ObsData%day=day
            ObsData%hour=hour
            ObsData%min=min
            ObsData%sec=int_sec
            ObsData%week=GPSweek
            ObsData%sow=GPSsec
            read(line, "(28X,I1,I3)") ObsData%Flag, ObsData%PRNS
            if ( allocated(ObsData%PRN) ) then
                deallocate(ObsData%PRN, ObsData%System, ObsData%C1,ObsData%C2, ObsData%P1, ObsData%P2, ObsData%P3, &
                        ObsData%L1, ObsData%L2, ObsData%L2C, ObsData%L3, ObsData%LLI1, ObsData%LLI2, ObsData%LLI2C, ObsData%LLI3, &
                        ObsData%D1, ObsData%D2, ObsData%D3)
            end if
            allocate(ObsData%PRN(ObsData%PRNS), ObsData%System(ObsData%PRNS), ObsData%C1(ObsData%PRNS), &
                        ObsData%C2(ObsData%PRNS),ObsData%P1(ObsData%PRNS), ObsData%P2(ObsData%PRNS),ObsData%P3(ObsData%PRNS), &
                        ObsData%L1(ObsData%PRNS), ObsData%L2(ObsData%PRNS), ObsData%L2C(ObsData%PRNS), ObsData%L3(ObsData%PRNS), &
                        ObsData%LLI1(ObsData%PRNS), ObsData%LLI2(ObsData%PRNS), ObsData%LLI2C(ObsData%PRNS), ObsData%LLI3(ObsData%PRNS),&
                        ObsData%D1(ObsData%PRNS), ObsData%D2(ObsData%PRNS), ObsData%D3(ObsData%PRNS) )
            ObsData%C1=0.0d0
            ObsData%C2=0.0d0
            ObsData%P1=0.0d0
            ObsData%P2=0.0d0
            ObsData%P3=0.0d0
            ObsData%L1=0.0d0
            ObsData%L2=0.0d0
            ObsData%L2C=0.0d0
            ObsData%L3=0.0d0
            ObsData%LLI1=0
            ObsData%LLI2=0
            ObsData%LLI2C=0
            ObsData%LLI3=0
            ObsData%D1=0.0d0
            ObsData%D2=0.0d0
            ObsData%D3=0.0d0

            line=line(33:68)
            if ( ObsData%PRNS>12) then
                read(unit=ObsID(sta),fmt="(32X,A)") line((len_trim(line)+1):)
                if ( ObsData%PRNS>24 ) then
                    read(unit=ObsID(sta),fmt="(32X,A)") line((len_trim(line)+1):)
                    if ( ObsData%PRNS>36 ) then
                        write(*,*) "Observation PRNS more than 36."
                        write(CoorID,*) "Stop in: "
                        write(CoorID,*) "Observation PRNS more than 36."
                        pause
                        stop
                    end if
                end if
            end if
        
            do i=1,ObsData%PRNS
                read(line(i*3-2:i*3),"(A1,I2)",end=200) ObsData%System(i), ObsData%PRN(i)
                if (ObsData%System(i)==" ")  ObsData%System(i)="G"
                ! Read the data
                read(unit=ObsID(sta), fmt="(A)",end=200) line2
                if (ObsHead(sta)%ObsTypes(1)>5) then
                    read(unit=ObsID(sta),fmt="(A)",end=200) line2(81:)
                    if (ObsHead(sta)%ObsTypes(1)>10) then
                        read(unit=ObsID(sta),fmt="(A)",end=200) line2(161:)
                        if (ObsHead(sta)%ObsTypes(1)>15) then
                            read(unit=ObsID(sta),fmt="(A)",end=200) line2(241:)
                            if (ObsHead(sta)%ObsTypes(1)>20) then
                                read(unit=ObsID(sta),fmt="(A)",end=200) line2(321:)
                                if (ObsHead(sta)%ObsTypes(1)>25) then
                                    write(*,*) "Obs Types more than 25!"
                                    write(CoorID, *) "Stop in: "
                                    write(CoorID, *) "Obs Types more than 25."
                                    pause
                                    stop
                                end if
                            end if
                        end if
                    end if
                end if
            
                do j=1,15   ! For all types of each satellite
                    k=ObsHead(sta)%ObsType(j,1)  ! Get the order of this type of observation
                    if (k==0) cycle  ! If this kind of observation not found
                    read(line2(16*k-15:16*k-2),"(F14.3)") val
                    read(line2(16*k-1:16*k-1),"(I1)") LLI
                
                    if (j==1) then
                        ObsData%C1(i)=val
                    else if (j==2) then
                        ObsData%P1(i)=val
                    else if (j==3) then
                        ObsData%P2(i)=val
                    else if (j==4) then
                        ObsData%P3(i)=val
                    else if (j==5) then
                        ObsData%L1(i)=val
                        ObsData%LLI1(i)=LLI  ! Loss of lock indicator
                    else if (j==6) then
                        ObsData%L2(i)=val
                        ObsData%LLI2(i)=LLI
                    else if (j==7) then
                        ObsData%L3(i)=val
                        ObsData%LLI3(i)=LLI
                    else if (j==8) then
                        ObsData%C2(i)=val
                    else if (j==9) then
                        if (val<LimSNR) then
                            ObsData%C1(i)=0.d0
                            ObsData%P1(i)=0.d0
                        end if
                    else if (j==10) then
                        if (val<LimSNR) ObsData%P2(i)=0.d0
                    else if (j==11) then
                        if (val<LimSNR) ObsData%P3(i)=0.d0
                    else if (j==12) then
                        ObsData%D1(i)=val
                    else if (j==13) then
                        ObsData%D2(i)=val
                    else if (j==14) then
                        ObsData%D3(i)=val
                    end if
                end do   !  do j=1,15
                
            end do  ! do i=1,ObsData%PRNS
        
            ! If next line not reach the specialized time
            if ((GPSweek-Obsweek)*604800.d0+GPSsec-Obssec<-0.01d0) then
!            if ((GPSweek<Obsweek) .or. (GPSsec<Obssec)) then
                cycle
            else
                exit
            end if
        end do   ! do while(.true.)
    
    else if(ObsHead(sta)%Version==3) then
        do while(.true.)
            ObsData%Flag=1   ! First set the flag false, in case that there is no data in this epoch
            read(unit=ObsID(sta),fmt="(A)",end=200) line
            if ( (index(line,"COMMENT") /=0) .or. (len_trim(line)==0) .or. (index(line,"SECONDS") /=0) ) then
                cycle
            end if
        
            read(line, "(31X,I1,I3)") Flag, lines
            if (Flag>1) then      ! Event flag
                do i=1,lines
                    read(unit=ObsID(sta), fmt="(A)",end=200) line
                    if (index(line(61:80),"ANTENNA: DELTA H/E/N") /=0 ) then  ! Event flag=3, eg. bjfs1660.10o
                        read(line,"(3F14.4)") ObsHead(sta)%Antenna(3), ObsHead(sta)%Antenna(2),ObsHead(sta)%Antenna(1)
                    end if
                end do
                cycle
            end if

            read(unit=line,fmt="(2X,I4,4(1X,I2),F11.7)") year, mon, day, hour, min, sec
            int_sec=real(nint(sec*10.d0))/10.d0   ! nint: nearest integer, for 1~10Hz data
            ObsData%Clk_Bias=sec-int_sec   ! station receiver clock correction
            ! int_sec=int_sec-14.d0  ! 精密星历时间为BDT，观测时间为GPST

            call UTC2GPST(year, mon, day, hour, min, int_sec, GPSweek, GPSsec)
            if ((GPSweek-Obsweek)*604800.d0+GPSsec-Obssec>0.01d0) then ! If no data in this epoch
            !if ((GPSweek>Obsweek) .or. (GPSsec>Obssec)) then ! If no data in this epoch
                backspace(ObsID(sta))
                exit
            end if

            if (year<80) then
                year=year+2000
            else if ((year>=80) .and. (year<100)) then
                year=year+1900
            end if
            ObsData%year=year
            ObsData%mon=mon
            ObsData%day=day
            ObsData%hour=hour
            ObsData%min=min
            ObsData%sec=int_sec
            ObsData%week=GPSweek
            ObsData%sow=GPSsec
            read(line, "(31X,I1,I3)") ObsData%Flag, ObsData%PRNS

            if ( allocated(ObsData%PRN) ) then
                deallocate(ObsData%PRN, ObsData%System, ObsData%C1, ObsData%C2, ObsData%P1, ObsData%P2,ObsData%P3, &
                        ObsData%L1, ObsData%L2, ObsData%L2C, ObsData%L3, ObsData%LLI1, ObsData%LLI2, ObsData%LLI2C, ObsData%LLI3, &
                        ObsData%D1, ObsData%D2, ObsData%D3)
            end if
            allocate(ObsData%PRN(ObsData%PRNS), ObsData%System(ObsData%PRNS), ObsData%C1(ObsData%PRNS), &
                        ObsData%C2(ObsData%PRNS),ObsData%P1(ObsData%PRNS), ObsData%P2(ObsData%PRNS),ObsData%P3(ObsData%PRNS), &
                        ObsData%L1(ObsData%PRNS), ObsData%L2(ObsData%PRNS), ObsData%L2C(ObsData%PRNS), ObsData%L3(ObsData%PRNS), &
                        ObsData%LLI1(ObsData%PRNS), ObsData%LLI2(ObsData%PRNS), ObsData%LLI2C(ObsData%PRNS), ObsData%LLI3(ObsData%PRNS), &
                        ObsData%D1(ObsData%PRNS), ObsData%D2(ObsData%PRNS), ObsData%D3(ObsData%PRNS) )
            ObsData%C1=0.0d0
            ObsData%C2=0.0d0
            ObsData%P1=0.0d0
            ObsData%P2=0.0d0
            ObsData%P3=0.0d0
            ObsData%L1=0.0d0
            ObsData%L2=0.0d0
            ObsData%L2C=0.0d0
            ObsData%L3=0.0d0
            ObsData%LLI1=0
            ObsData%LLI2=0
            ObsData%LLI2C=0
            ObsData%LLI3=0
            ObsData%D1=0.0d0
            ObsData%D2=0.0d0
            ObsData%D3=0.0d0

            do i=1,ObsData%PRNS
                read(unit=ObsID(sta), fmt="(A1,I2,A)",end=200) System,PRN, line2
                if ( (index(line,"COMMENT") /=0) .or. (len_trim(line2)==0) .or. (index(line,"SECONDS") /=0) ) then
                    cycle
                end if
                ObsData%System(i)=system
                ObsData%PRN(i)=PRN
                if (ObsData%System(i)==" ")  ObsData%System(i)="G"
                if (ObsData%System(i)=="G") then
                    sys=1
                elseif (ObsData%System(i)=="R") then
                    sys=2
                else if (ObsData%System(i)=="C") then
                    sys=3
                else if (ObsData%System(i)=="E") then
                    sys=4
                else if (ObsData%System(i)=="J") then
                    sys=5
                else if (ObsData%System(i)=="I") then
                    sys=6
                else
                    cycle
                end if

                do j=1,15   ! For all types of each satellite
                    k=ObsHead(sta)%ObsType(j,sys)  ! Get the order of this type of observation
                    if (k==0) cycle  ! If this kind of observation not found
                    read(line2(16*k-15:16*k-2),"(F14.3)") val
                    read(line2(16*k-1:16*k-1),"(I1)") LLI
                
                    if (j==1) then
                        ObsData%C1(i)=val
                    else if (j==2) then
                        ObsData%P1(i)=val
                    else if (j==3) then
                        ObsData%P2(i)=val
                    else if (j==4) then
                        ObsData%P3(i)=val
                    else if (j==5) then
                        ObsData%L1(i)=val
                        ObsData%LLI1(i)=LLI
                    else if (j==6) then
                        ObsData%L2(i)=val
                        ObsData%LLI2(i)=LLI
                    else if (j==7) then
                        ObsData%L3(i)=val
                        ObsData%LLI3(i)=LLI
                    else if (j==8) then
                        ObsData%C2(i)=val
                    else if (j==9) then
                        if (val<LimSNR) then
                            ObsData%C1(i)=0.d0
                            ObsData%P1(i)=0.d0
                        end if
                    else if (j==10) then
                        if (val<LimSNR) ObsData%P2(i)=0.d0
                    else if (j==11) then
                        if (val<LimSNR) ObsData%P3(i)=0.d0
                    else if (j==12) then
                        ObsData%D1(i)=val
                    else if (j==13) then
                        ObsData%D2(i)=val
                    else if (j==14) then
                        ObsData%D3(i)=val
                    else if (j==15) then
                        ObsData%L2C(i)=val
                        ObsData%LLI2C(i)=LLI
                    end if
                end do   !  do j=1,5

            end do   ! do i=1,ObsData%PRNS

            ! If next line not reach the specialized time
            if ((GPSweek-Obsweek)*604800.d0+GPSsec-Obssec<-0.01d0) then
!            if ((GPSweek<Obsweek) .or. (GPSsec<Obssec)) then
                cycle
            else
                exit
            end if
        end do
    end if
    200 return ! When reach the end of the file
    
end subroutine