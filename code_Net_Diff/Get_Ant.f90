! ========== Get Antenna Information ========
!
! PURPOSE:
!      Read the antenna information file 'igs08_****.atx' 
!  and get the PCO and PCVof satellite and receiver antenna
!  at a given MJD time.
!
!  REFERENCES:
!        http://www.epncb.oma.be/ftp/station/general/antex14.txt
!        http://www.navipedia.net/index.php/Antenna_Phase_Centre
!        http://www.navipedia.net/index.php/Satellite_Antenna_Phase_Centre
!        http://www.navipedia.net/index.php/Receiver_Antenna_Phase_Centre
!
! INPUTS:
!      MJD             Modified Julian Day
!      Ant_Type     Antenna type
!
! OUTPUT(as a public module data):
!      Ant             A type data, See MOD.F90
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!   ============ End of Header ==========

subroutine Get_Ant(MJD)
use MOD_Ant
use MOD_constant
use MOD_FileID      ! AntID
use MOD_FileDir    ! AntFile
use MOD_STA
implicit none
    logical :: alive, IsRecAnt=.false.
    real(8) :: MJD, ValidMJD
!    character(20) :: PCOFile
    character(20) :: Sta_Ant
    ! Local variables
    character(400) :: line
    character(30) :: keyword
    character(20) :: Ant_Type 
    character :: System
    integer :: PRN, i=0, j, k, l,Flag
    real(8) :: dazi, zen1, zen2, dzen
    integer :: nazi, nzen
    integer :: year, mon, day, hour, min
    real(8) :: sec, temp
    
    Flag=0
    inquire(file=AntFile,exist=alive)
    if (.not. alive) then
        write(*,*) "Antenna file: """//trim(AntFile)//""" doesn't exist!"
        pause
        stop
    end if
    AntID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=AntID, file=AntFile)
    
    do while(.true.)
        read(unit=AntID,fmt='(A)',end=200) line
        keyword=line(61:80)
        if (index(keyword,"END OF HEADER") /= 0)  exit
    end do
    
    do while(.true.)
        IsRecAnt=.false.
        read(unit=AntID,fmt='(A)',end=200) line
        keyword=line(61:80)
        if  (index(keyword,"START OF ANTENNA") /= 0) then
            do while(.true.)
                read(unit=AntID,fmt='(A)',end=200) line
                keyword=line(61:80)
                if  (index(keyword,"TYPE / SERIAL NO") /= 0) then
                    read(line,'(A20,A1,I2)') Ant_Type, System, PRN
                    if (System=='G') then
                        if (PRN>GNum) exit
                        Ant(PRN)%PRN=PRN
                        Ant(PRN)%Ant_Type=Ant_Type
                    else if (System=='R') then
                        if (PRN>RNum) exit
                        PRN=PRN+GNum
                        Ant(PRN)%PRN=PRN
                        Ant(PRN)%Ant_Type=Ant_Type
                    else if (System=='C') then
                        if (PRN>CNum) exit
                        PRN=PRN+GNum+RNum
                        Ant(PRN)%PRN=PRN
                        Ant(PRN)%Ant_Type=Ant_Type
                    else if (System=='E') then
                        if (PRN>NumE) exit
                        PRN=PRN+GNum+RNum+CNum
                        Ant(PRN)%PRN=PRN
                        Ant(PRN)%Ant_Type=Ant_Type
                    else if (System=='J') then
                        if (PRN>JNum) exit
                        PRN=PRN+GNum+RNum+CNum+NumE
                        Ant(PRN)%PRN=PRN
                        Ant(PRN)%Ant_Type=Ant_Type
                    else if (System=='I') then
                        if (PRN>INum) exit
                        PRN=PRN+GNum+RNum+CNum+NumE+JNum
                        Ant(PRN)%PRN=PRN
                        Ant(PRN)%Ant_Type=Ant_Type
                    else
                        do i=1, STA%Num
                            if ((STA%STA(i)%Ant_Type==Ant_Type) .and. (Ant(GNum+RNum+CNum+NumE+JNum+i)%PRN==0)) then
                                PRN=GNum+RNum+CNum+NumE+JNum+i
                                Ant(PRN)%PRN=PRN
                                Ant(PRN)%Ant_Type=Ant_Type
                                Flag=Flag+1
                                IsRecAnt=.true.
                                exit
                            end if
                        end do
                        if (IsRecAnt==.false.) exit
                    end if
                else if (index(keyword,"DAZI") /= 0) then   ! Azimuth
                    read(line,'(2X,F6.1)') dazi
                    if (dazi==0.d0) then ! non-azimuth-dependent
                        nazi=0
                    else
                        nazi=nint(360.d0/dazi)+1
                    end if
                else if (index(keyword,"ZEN1 / ZEN2 / DZEN") /= 0) then   ! Zenith
                    read(line,'(2X,3F6.1)') zen1, zen2, dzen
                    nzen=nint((zen2-zen1)/dzen)+1
                    
                else if (index(keyword,"# OF FREQUENCIES") /= 0) then
                    read(line,'(I6)') j    ! Get number of frequency
                    if (j<2) j=2   ! If only one frequency, then set another as zero
                    if (System=='E' .and. j<5) j=5 ! For Galileo
                    k=0
                else if (index(keyword,"VALID FROM") /= 0) then
                    read(line,'(5I6,F13.7)') year,mon, day, hour, min, sec
                    call UTC2MJD(year, mon, day, hour, min, sec, ValidMJD)
                    if (MJD<ValidMJD) exit   !  Exit a type of antenna
                else if (index(keyword,"VALID UNTIL") /= 0) then
                    read(line,'(5I6,F13.7)') year,mon, day, hour, min, sec
                    call UTC2MJD(year, mon, day, hour, min, sec, ValidMJD)
                    if (MJD>ValidMJD) exit   !  Exit a type of antenna
                else if (index(keyword,"START OF FREQUENCY") /= 0) then
                    if (.not. allocated(Ant(PRN)%Freq)) then
                         allocate(Ant(PRN)%PCO(3,j))
                         Ant(PRN)%PCO=0.d0
                         allocate(Ant(PRN)%Freq(j))
                    end if
                    if (.not. allocated(Ant(PRN)%PCV)) then
                        allocate(Ant(PRN)%PCV(nzen,nazi+1,j))
                    end if
                    Ant(PRN)%nzen=nzen
                    Ant(PRN)%nazi=nazi
                    Ant(PRN)%dzen=dzen
                    Ant(PRN)%dazi=dazi
                    if (.not. allocated(Ant(PRN)%Zenith)) then
                        allocate(Ant(PRN)%Zenith(nzen))
                        allocate(Ant(PRN)%Azimuth(nazi))
                        do i=1,nzen ! zenith
                            Ant(PRN)%Zenith(i)=zen1+dzen*(i-1)
                        end do
                        do i=1,nazi ! azimuth
                            Ant(PRN)%Azimuth(i)=0.d0+dazi*(i-1)
                        end do
                    end if
                    if (System=='E') then
                        if (index(line(1:7),"E01") /= 0) then
                            k=1  ! E1
                        elseif (index(line(1:7),"E05") /= 0) then
                            k=2  ! E5a
                        elseif (index(line(1:7),"E07") /= 0) then
                            k=3  ! E5b
                        elseif (index(line(1:7),"E08") /= 0) then
                            k=4  ! E5
                        elseif (index(line(1:7),"E06") /= 0) then
                            k=5  ! E6
                        end if
                    else
                        k=k+1
                    end if
                    read(line,'(A6)') Ant(PRN)%Freq(k)
                    do while(.true.)
                        read(unit=AntID,fmt='(A)',end=200) line
                        keyword=line(61:80)
                        if (index(keyword,"NORTH / EAST / UP") /= 0) then
                            read(line,'(3F10.2)') Ant(PRN)%PCO(:,k)
                            Ant(PRN)%PCO(:,k)=Ant(PRN)%PCO(:,k)/1000.d0
                        else if (line(4:8)=="NOAZI") then    
                            read(line(9:),*) Ant(PRN)%PCV(1:nzen,1,k)
                            do l=1,nazi
                                read(AntID,*) temp, Ant(PRN)%PCV(1:nzen,l+1,k)
                            end do
                            Ant(PRN)%PCV(:,:,k)=Ant(PRN)%PCV(:,:,k)*1.d-3
                        else if (index(keyword,"END OF FREQUENCY") /= 0) then
                            exit
                        end if
                    end do
                end if   ! if  (index(keyword,"TYPE / SERIAL NO") /= 0)
                if (index(keyword,"END OF ANTENNA") /= 0) exit
            end do  ! do while(.true.)
        end if  ! if  (index(keyword,"START OF ANTENNA") /= 0)
        if (Flag==STA%Num) then
            exit  ! Exit read atx file
        elseif (IsRecAnt) then
            do while (.true.)
                backspace(AntID)   ! 给相同天线的不同测站的PCO/PCV赋值
                read(unit=AntID,fmt='(A)',end=200) line
                backspace(AntID)
                keyword=line(61:80)
                if  (index(keyword,"START OF ANTENNA") /= 0) exit
            end do
        end if
    end do  ! do while(.true.)
    
    200 close(AntID)
    do i=1,STA%Num
        if (Ant(GNum+RNum+CNum+NumE+JNum+INum+i)%PRN==0) then
            write(*,*) "Can't find the antenna information of the receiver at station "//STA%STA(i)%Name//", it will be set as 0."
            Ant(GNum+RNum+CNum+NumE+JNum+INum+i)%PRN=GNum+RNum+CNum+NumE+JNum+INum+i
            allocate(Ant(GNum+RNum+CNum+NumE+JNum+INum+i)%PCO(3,2))
            Ant(GNum+RNum+CNum+NumE+JNum+INum+i)%PCO=0.d0
            allocate(Ant(GNum+RNum+CNum+NumE+JNum+INum+i)%Freq(2))
!            pause
        end if
    end do
    return
end subroutine