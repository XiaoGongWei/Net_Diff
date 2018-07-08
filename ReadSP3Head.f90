! ===============  ReadSP3Data  =============
!
! PURPOSE:
!     Read the Header of SP3 File
!
! INPUTS:
!      n              read n epochs of sp3 data
!
! OUTPUT:
!      SP3Data   the added sp3 data, in SP3Data module
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!===============  End of Header  ==============

subroutine ReadSP3Head( SP3File)
use MOD_FileID
use MOD_SP3Head  ! SP3Head is a common and type variable
use MOD_VAR
implicit none
    character(len=200) :: SP3File
    logical :: alive
    integer :: error
    integer ::  i, j
    integer :: PRN
    character(len=400) :: line
    character(len=100) :: line2
    character(3) :: temp
        
    i=0
    j=1
    inquire(file=SP3File,exist=alive)
    if (.not. alive) then
        write(*,*) "SP3 flie: """//trim(SP3File)//""" doesn't exist!"
        pause
        stop
    end if
   

    SP3ID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=SP3ID, file=SP3File)
    
    
    ! Start to read SP3 file header
    ! Line 1
    read(unit=SP3ID,fmt=*,iostat=error)
    i=i+1
    ! Line 2
    read(SP3ID,fmt="(3X,I4,F15.8)",iostat=error) SP3Head%GPSWeek, SP3Head%GPSSec
    i=i+1
    ! Line 3, read PRNS and PRN
    read(unit=SP3ID,fmt="(4X,I2,3X,A51)",iostat=error) SP3Head%PRNS, line
    i=i+1
    ! Check satellite number, if > 17,34,51, read next line inorder. Till now, no more than 68.
    if (SP3Head%PRNS>17) then
        ! Line 4
        read(unit=SP3ID,fmt="(9X,A51)",iostat=error) line2
        line=trim(line)//trim(line2)
        i=i+1
        if (SP3Head%PRNS>34) then
            ! Line 5
            read(unit=SP3ID,fmt="(9X,A51)",iostat=error) line2
            line=trim(line)//trim(line2)
            i=i+1
            if (SP3Head%PRNS>51) then
                ! Line 6
                read(unit=SP3ID,fmt="(9X,A51)",iostat=error) line2
                line=trim(line)//trim(line2)
                i=i+1
                if (SP3Head%PRNS>68) then
                    ! Line 7
                    read(unit=SP3ID,fmt="(9X,A51)",iostat=error) line2
                    line=trim(line)//trim(line2)
                    i=i+1
                    if (SP3Head%PRNS>85) then
                        ! Line 8
                        read(unit=SP3ID,fmt="(9X,A51)",iostat=error) line2
                        line=trim(line)//trim(line2)
                        i=i+1
                    end if
                end if
            end if
        end if
    end if
    
    ! Set system and prn for each satellite
    allocate(SP3Head%SYSTEM(SP3Head%PRNS))
    allocate(SP3Head%PRN(SP3Head%PRNS))
    do j=1,SP3Head%PRNS
        read(line((j*3-2):(j*3)), fmt="(A1I2)") SP3Head%SYSTEM(j), SP3Head%PRN(j)
    end do
    
    line=''
    do while(.true.)
        read(SP3ID,fmt="(A80)",iostat=error) line2
        i=i+1
!        if (i==22) exit
        if (error /=0) exit
        if (index(line2,"  cc ") /= 0)   then  ! adapted to sp3_d
!        if (i==13) then
            read(line2,'(9X,A3)') temp
            if (temp=='GPS') then
                SP3Time=0.d0
            elseif ((temp=='BDS') .or. (temp=='BDT'))then
                SP3Time= -14.d0
            else 
                write(*,*) "Can't indentify string "//temp//" at line 13 in SP3 file"//SP3File
                pause
                stop
            end if
        end if
        if (index(line2,"++") /= 0)   then
            line=trim(line)//trim(line2(10:))
        elseif (index(line2(1:2),"* ") /= 0)   then
            backspace(SP3ID)
            exit
        end if
    end do

    do j=1,SP3Head%PRNS
        PRN=SP3Head%PRN(j)
        if (SP3Head%SYSTEM(j)=='R') then
            PRN=PRN+GNum
        else if (SP3Head%SYSTEM(j)=='C') then
            PRN=PRN+GNum+RNum
        else if (SP3Head%SYSTEM(j)=='E') then
            PRN=PRN+GNum+RNum+CNum
        else if (SP3Head%SYSTEM(j)=='J') then
            PRN=PRN+GNum+RNum+CNum+NumE
        else if (SP3Head%SYSTEM(j)=='I') then
            PRN=PRN+GNum+RNum+CNum+NumE+JNum
        elseif (SP3Head%SYSTEM(j)=="G") then
            prn=prn
        else
            cycle
        end if
        read(line((j*3-1):(j*3)), fmt="(I2)") SP3Head%OA(PRN)
    end do

    return
end subroutine