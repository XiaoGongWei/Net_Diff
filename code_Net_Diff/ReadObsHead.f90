! Read Observation File Header
! ================  ReadObsData ==============
!
! PURPOSE:
!     Read Observation File Header from **.**o file
!
! NOTICE:
!     1) The time system should be consistent with the
!     orbit and clock.
!     2) Both RINEX3.02 and 2.01 format can be handled here
!         For more information about RINEX file, refer to:
!            http://spot.colorado.edu/~kristine/rinex210.txt
!            ftp://igs.org/pub/data/format/rinex302.pdf
!            Steigenberger P, Montenbruck O, Hessel U. Performance 
!         Evaluation of the Early CNAV Navigation Message[J]. Navigation, 
!         2015, 62(3):219-228.
!     3) For RINEX 3.01 and 3.02, the '# / TYPES OF OBSERV'
!          is not sure for there are many types in every type.
!          This to be further improved.
!
! INPUTS:
!      ObsFile        observation file name
!      kth               order of multi-stations
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================  End of Header  ===============

subroutine ReadObsHead(ObsFile,kth)
use MOD_FileID
use MOD_ObsHead
use MOD_VAR
use MOD_STA
implicit none
    character(len=200) ObsFile
    logical :: alive
    integer :: i, j, kth
    character(len=300) :: line
    character(len=81) :: line2
    character(len=30) keyword
    character(len=2) :: TempType
    character(len=3) :: TempType3
    integer :: year, mon,day, hour, min
    real(8) :: sec, MJD1, MJD2
    integer(1) :: ObsTypes

    ! inquire that if the obs file exists
    inquire(file=ObsFile,exist=alive)
    if (.not. alive) then
        write(*,*) "Observation flie: """//trim(ObsFile)//""" doesn't exist!"
        pause
        stop
    end if
    
    ObsID(kth)=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=ObsID(kth), file=ObsFile,action="read")
    
    Read_File_Loop: do while(.true.)
        read(unit=ObsID(kth),fmt="(A)") line
        keyword=line(61:80)
        
        if (index(keyword,"RINEX VERSION / TYPE") /= 0) then
            read(line,fmt="(5X,I1)") ObsHead(kth)%Version
            if (ObsHead(kth)%Version==2) then
                allocate(ObsHead(kth)%ObsTypes(1))
                ObsHead(kth)%ObsTypes=0
                allocate(ObsHead(kth)%ObsType(15,1))
                ObsHead(kth)%ObsType=0
            else if (ObsHead(kth)%Version==3) then
                allocate(ObsHead(kth)%ObsTypes(6))
                ObsHead(kth)%ObsTypes=0
                allocate(ObsHead(kth)%ObsType(15,6))
                ObsHead(kth)%ObsType=0
            end if
        elseif (index(keyword,"# / TYPES OF OBSERV") /= 0) then
            read(line,*) ObsHead(kth)%ObsTypes(1)
            if (ObsHead(kth)%ObsTypes(1) > 9) then
                read(unit=ObsID(kth),fmt="(A)") line2
                line=trim(line(1:60))//trim(line2(7:60))
                if (ObsHead(kth)%ObsTypes(1) > 18) then
                    read(unit=ObsID(kth),fmt="(A)") line2
                    line=trim(line)//trim(line2(7:60))
                    if (ObsHead(kth)%ObsTypes(1)> 27) then
                        write(*,*) "ObsTypes more than 27."
                        write(CoorID,*) "Stop in:"
                        write(CoorID,*) " Obs types more than 27."
                        pause
                        stop
                    end if
                end if
            end if
            Find_Type: do j=1,ObsHead(kth)%ObsTypes(1)
                TempType=line((j*6+5):(j*6+6))
                if (TempType=="C1") then
                    ObsHead(kth)%ObsType(1,1)=j
                else if (TempType=="P1") then
                    ObsHead(kth)%ObsType(2,1)=j
                else if (TempType== "P2") then
                    ObsHead(kth)%ObsType(3,1)=j
                else if (TempType== "P3") then
                    ObsHead(kth)%ObsType(4,1)=j
                else if (TempType=="L1") then
                    ObsHead(kth)%ObsType(5,1)=j
                else if (TempType== "L2") then
                    ObsHead(kth)%ObsType(6,1)=j
                else if (TempType== "L3") then
                    ObsHead(kth)%ObsType(7,1)=j
                else if (TempType== "C2") then
                    ObsHead(kth)%ObsType(8,1)=j
                else if (TempType== "S1") then
                    ObsHead(kth)%ObsType(9,1)=j
                else if (TempType== "S2") then
                    ObsHead(kth)%ObsType(10,1)=j
                else if (TempType== "S3") then
                    ObsHead(kth)%ObsType(11,1)=j
                else if (TempType== "D1") then
                    ObsHead(kth)%ObsType(12,1)=j
                else if (TempType== "D2") then
                    ObsHead(kth)%ObsType(13,1)=j
                else if (TempType== "D3") then
                    ObsHead(kth)%ObsType(14,1)=j
                end if
            end do Find_Type
        elseif (index(keyword,"SYS / # / OBS TYPES") /= 0) then
            read(line(5:7),*) ObsTypes
            if (ObsTypes> 13) then
                read(unit=ObsID(kth),fmt="(A)") line2
                line=trim(line(1:60))//trim(line2(7:60))
                if (ObsTypes > 26) then
                    read(unit=ObsID(kth),fmt="(A)") line2
                    line=trim(line)//trim(line2(7:60))
                    if (ObsTypes > 39) then
                        read(unit=ObsID(kth),fmt="(A)") line2
                        line=trim(line)//trim(line2(7:60))
                        if (ObsTypes > 52) then
                            write(*,*) "ObsTypes more than 52."
                            write(CoorID,*) "Stop in:"
                            write(CoorID,*) " Obs types more than 52."
                            pause
                            stop
                        end if
                    end if
                end if
            end if
            if (line(1:1)=="G") then
                i=1
            elseif (line(1:1)=="R") then
                i=2
            elseif (line(1:1)=="C") then
                i=3
            elseif (line(1:1)=="E") then
                i=4
            elseif (line(1:1)=="J") then
                i=5
            elseif (line(1:1)=="I") then
                i=6
            else
                cycle
            end if
            ObsHead(kth)%ObsTypes(i)=ObsTypes
            Find_Type_2: do j=1,ObsHead(kth)%ObsTypes(i)
                TempType3=line((j*4+4):(j*4+6))
                if (i==1) then  ! GPS
                    if (TempType3=="C1C") then  ! L1C/A, C/A code
                        ObsHead(kth)%ObsType(1,i)=j
                    else if (TempType3=="C1W") then  ! P1(Y), Z-tracking and similar
                        ObsHead(kth)%ObsType(2,i)=j
                    else if (TempType3=="C2W") then  ! P2(Y), Z-tracking and similar
                        ObsHead(kth)%ObsType(3,i)=j
                    else if (TempType3=="C5Q" .or. TempType3=="C5X" ) then   ! P3
                        ObsHead(kth)%ObsType(4,i)=j
                    else if (TempType3=="L1C") then   ! L1
                        ObsHead(kth)%ObsType(5,i)=j
                    else if (TempType3=="L2W") then  ! L2
                        ObsHead(kth)%ObsType(6,i)=j
                    else if (TempType3=="L5Q" .or. TempType3=="L5X" ) then   ! L3
                        ObsHead(kth)%ObsType(7,i)=j
                    else if (TempType3=="C2S") then   ! L2C(M), L2C moderate length code
                        ObsHead(kth)%ObsType(8,i)=j
                    else if (TempType3=="C2L") then   ! L2C(L), L2C longer length code
                        ObsHead(kth)%ObsType(8,i)=j
                    else if (TempType3=="C2X") then   ! L2C(M+L), combined tracking of C2S and C2L
                        ObsHead(kth)%ObsType(8,i)=j
                    else if (TempType3=="S1C") then   ! S1
                        ObsHead(kth)%ObsType(9,i)=j
                    else if (TempType3=="S2W") then   ! S2
                        ObsHead(kth)%ObsType(10,i)=j
                    else if (TempType3=="S5Q") then   ! S3
                        ObsHead(kth)%ObsType(11,i)=j
                    else if (TempType3=="D1C") then   ! D1
                        ObsHead(kth)%ObsType(12,i)=j
                    else if (TempType3=="D2W") then   ! D2
                        ObsHead(kth)%ObsType(13,i)=j
                    else if (TempType3=="D5Q") then   ! D3
                        ObsHead(kth)%ObsType(14,i)=j
                    else if (TempType3=="L2S") then   ! L2C(M), L2C moderate length code
                        ObsHead(kth)%ObsType(15,i)=j
                    else if (TempType3=="L2L") then   ! L2C(L), L2C longer length code
                        ObsHead(kth)%ObsType(15,i)=j
                    else if (TempType3=="L2X") then   ! L2C(M+L), combined tracking of C2S and C2L
                        ObsHead(kth)%ObsType(15,i)=j
                    end if
                elseif (i==2) then   ! GLONASS
                    if (TempType3=="C1C") then
                        ObsHead(kth)%ObsType(1,i)=j
                    else if (TempType3=="C1P") then
                        ObsHead(kth)%ObsType(2,i)=j
                    else if (TempType3=="C2P") then
                        ObsHead(kth)%ObsType(3,i)=j
                    else if (TempType3=="C3I") then
                        ObsHead(kth)%ObsType(4,i)=j
                    else if (TempType3=="L1C") then
                        ObsHead(kth)%ObsType(5,i)=j
                    else if (TempType3=="L2C") then
                        ObsHead(kth)%ObsType(6,i)=j
                    else if (TempType3=="L3I") then
                        ObsHead(kth)%ObsType(7,i)=j
                    else if (TempType3=="C2C") then
                        ObsHead(kth)%ObsType(8,i)=j
                    else if ( (TempType3=="S1C") .or. (TempType3=="S1P") ) then
                        ObsHead(kth)%ObsType(9,i)=j
                    else if ( (TempType3=="S2C") .or. (TempType3=="S2P") ) then
                        ObsHead(kth)%ObsType(10,i)=j
                    else if (TempType3=="S3I") then
                        ObsHead(kth)%ObsType(11,i)=j
                    else if ( (TempType3=="D1C") .or. (TempType3=="D1P") ) then   ! D1
                        ObsHead(kth)%ObsType(12,i)=j
                    else if ( (TempType3=="D2C") .or. (TempType3=="D2P") ) then   ! D2
                        ObsHead(kth)%ObsType(13,i)=j
                    else if ( (TempType3=="D3I") .or. (TempType3=="D3X") ) then   ! D3
                        ObsHead(kth)%ObsType(14,i)=j
                    end if
                elseif (i==3) then   ! COMPASS
                    if (TempType3=="Q2") then
                        ObsHead(kth)%ObsType(1,i)=j
                    else if ( (TempType3=="C1I") .or. (TempType3=="C2I") .or. (TempType3=="C1X") .or. (TempType3=="C2X") ) then
                        ObsHead(kth)%ObsType(2,i)=j    ! C2I is for RINEX 3.01/3.03 and C1I is for RINEX 3.02
                    else if ( (TempType3=="C7I") .or. (TempType3=="C7X") ) then
                        ObsHead(kth)%ObsType(3,i)=j
                    else if ( (TempType3=="C6I") .or. (TempType3=="C6X") ) then
                        ObsHead(kth)%ObsType(4,i)=j
                    else if ( (TempType3=="L1I") .or. (TempType3=="L2I") .or. (TempType3=="L1X") .or. (TempType3=="L2X")) then
                        ObsHead(kth)%ObsType(5,i)=j  ! L2I is for RINEX 3.01/3.03 and L1I is for RINEX 3.02
                    else if ( (TempType3=="L7I") .or. (TempType3=="L7X") ) then
                        ObsHead(kth)%ObsType(6,i)=j
                    else if ( (TempType3=="L6I") .or. (TempType3=="L6X") ) then
                        ObsHead(kth)%ObsType(7,i)=j
                    else if ( (TempType3=="S1I") .or. (TempType3=="S2I") .or. (TempType3=="S1X") .or. (TempType3=="S2X") ) then
                        ObsHead(kth)%ObsType(9,i)=j
                    else if ( (TempType3=="S7I")  .or. (TempType3=="S7X") ) then
                        ObsHead(kth)%ObsType(10,i)=j
                    else if ( (TempType3=="S6I") .or. (TempType3=="S6X") ) then
                        ObsHead(kth)%ObsType(11,i)=j
                    else if ( (TempType3=="D1I") .or. (TempType3=="D1X") .or. (TempType3=="D2I") .or. (TempType3=="D2X") ) then   ! D1
                        ObsHead(kth)%ObsType(12,i)=j
                    else if ( (TempType3=="D7I") .or. (TempType3=="D7X") ) then   ! D2
                        ObsHead(kth)%ObsType(13,i)=j
                    else if ( (TempType3=="D6I") .or. (TempType3=="D6X") ) then   ! D3
                        ObsHead(kth)%ObsType(14,i)=j
                    end if
                elseif (i==4) then   ! GALILEO
                    if (TempType3=="Q2") then
                        ObsHead(kth)%ObsType(1,i)=j
                    else if ( (TempType3=="C1X") .or. (TempType3=="C1C") ) then
                        ObsHead(kth)%ObsType(2,i)=j
                    else if ( (TempType3=="C5X") .or. (TempType3=="C5Q") ) then
                        ObsHead(kth)%ObsType(3,i)=j
                    else if ( (TempType3=="C7X") .or. (TempType3=="C7Q") ) then
                        ObsHead(kth)%ObsType(4,i)=j
                    else if ( (TempType3=="L1X") .or. (TempType3=="L1C") ) then
                        ObsHead(kth)%ObsType(5,i)=j
                    else if ( (TempType3=="L5X") .or. (TempType3=="L5Q") ) then
                        ObsHead(kth)%ObsType(6,i)=j
                    else if ( (TempType3=="L7X") .or. (TempType3=="L7Q") ) then
                        ObsHead(kth)%ObsType(7,i)=j
                    else if ( (TempType3=="S1X") .or. (TempType3=="S1C") ) then
                        ObsHead(kth)%ObsType(9,i)=j
                    else if ( (TempType3=="S5X") .or. (TempType3=="S5Q") ) then
                        ObsHead(kth)%ObsType(10,i)=j
                    else if ( (TempType3=="S7X") .or. (TempType3=="S7Q") ) then
                        ObsHead(kth)%ObsType(11,i)=j
                    else if ( (TempType3=="D1C") .or. (TempType3=="D1X") ) then   ! D1
                        ObsHead(kth)%ObsType(12,i)=j
                    else if ( (TempType3=="D5Q") .or. (TempType3=="D5X") ) then   ! D2
                        ObsHead(kth)%ObsType(13,i)=j
                    else if ( (TempType3=="D7Q") .or. (TempType3=="D7X") ) then   ! D3
                        ObsHead(kth)%ObsType(14,i)=j
                    end if
                elseif (i==5) then   ! QZSS
                    if (TempType3=="C1C") then
                        ObsHead(kth)%ObsType(1,i)=j
                    else if ( (TempType3=="C1X")  .or. (TempType3=="C1Z") ) then
                        ObsHead(kth)%ObsType(2,i)=j
                    else if ( (TempType3=="C2X") ) then
                        ObsHead(kth)%ObsType(3,i)=j
                    else if ( (TempType3=="C5X") .or. (TempType3=="C5Q") .or. (TempType3=="C5I") ) then
                        ObsHead(kth)%ObsType(4,i)=j
                    else if (TempType3=="L1C") then  ! L1X L1Z it is the new signal
                        ObsHead(kth)%ObsType(5,i)=j
                    else if ( (TempType3=="L2X") ) then
                        ObsHead(kth)%ObsType(6,i)=j
                    else if ( (TempType3=="L5X") .or. (TempType3=="L5Q") .or. (TempType3=="L5I") ) then
                        ObsHead(kth)%ObsType(7,i)=j
                    else if ( (TempType3=="C2L") .or. (TempType3=="C2S") ) then
                        ObsHead(kth)%ObsType(8,i)=j
                    else if ( (TempType3=="S1X") .or. (TempType3=="S1C")  .or. (TempType3=="S1Z") ) then
                        ObsHead(kth)%ObsType(9,i)=j
                    else if ( (TempType3=="S2X") .or. (TempType3=="S2L")  .or. (TempType3=="S2S") ) then
                        ObsHead(kth)%ObsType(10,i)=j
                    else if ( (TempType3=="S5X") .or. (TempType3=="S5Q")  .or. (TempType3=="S5I") ) then
                        ObsHead(kth)%ObsType(11,i)=j
                    else if ( (TempType3=="D1X") .or. (TempType3=="D1C") ) then   ! D1
                        ObsHead(kth)%ObsType(12,i)=j
                    else if ( (TempType3=="D2X") .or. (TempType3=="D2L") ) then   ! D2
                        ObsHead(kth)%ObsType(13,i)=j
                    else if ( (TempType3=="D5X") .or. (TempType3=="D5Q") ) then   ! D3
                        ObsHead(kth)%ObsType(14,i)=j
                    else if ( (TempType3=="L2L") .or. (TempType3=="L2S") ) then
                        ObsHead(kth)%ObsType(15,i)=j
                    end if
                elseif (i==6) then   ! IRNSS
                    if (TempType3=="C5A") then
                        ObsHead(kth)%ObsType(2,i)=j
                    else if (TempType3=="L5A") then
                        ObsHead(kth)%ObsType(5,i)=j
                    else if (TempType3=="S5A") then
                        ObsHead(kth)%ObsType(9,i)=j
                    end if
                else
                    write(*,*) 'not defined'
                    pause
                end if
            end do Find_Type_2
        else if (index(keyword,"APPROX POSITION XYZ") /= 0) then
            read(line,*)  ObsHead(kth)%AppCoor(1:3)
        else if (index(keyword,"ANT # / TYPE") /= 0) then
            read(line,'(20X,A20)') ObsHead(kth)%Ant_Type
        else if (index(keyword,"ANTENNA: DELTA H/E/N") /= 0) then
            read(line, *)  ObsHead(kth)%Antenna(3), ObsHead(kth)%Antenna(2), ObsHead(kth)%Antenna(1)
            STA%STA(kth)%NEU=ObsHead(kth)%Antenna
            !STA%Coor(1:3,i)=STA%Coor(1:3,i)+MATMUL(ObsHead(kth)%Antenna, STA%Rota(i)%Rotation)
        else if (index(keyword,"INTERVAL") /= 0) then
            read(line,*)  ObsHead(kth)%Interval
            if (kth>STA%FixNum .and. Interval<ObsHead(kth)%Interval) Interval=ObsHead(kth)%Interval
        else if (index(keyword,"TIME OF FIRST OBS") /= 0) then
            read(line,'(5I6,F13.7,5X,A3)') year, mon, day, hour, min, sec, TimeSys
            if (TimeSys=="GLO") then
                write(*,*) "System time is in GLONASS time in Observation Header."
                pause
                stop
            end if
            call UTC2GPST(year, mon, day, hour, min, sec, ObsHead(kth)%GPSweek, ObsHead(kth)%GPSsec)
        else if (index(keyword,"END OF HEADER") /= 0) then
            exit
        end if
    end do Read_File_Loop
    return
end subroutine