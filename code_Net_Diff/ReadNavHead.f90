! ===============  ReadNavHead  =============
!
! PURPOSE:
!    Read the Header of GPS navigation ephemeris file
!
! INPUT:
!    NavFile            navigation ephemeris file
!  
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!  ===============  End of Header  =============

subroutine ReadNavHead(NavFile)
use MOD_FileID
use MOD_NavHead
implicit none
    character(200) :: NavFile
    logical :: alive
    integer :: error
    character(100) :: line
    character(20) :: keyword
    
    inquire(file=NavFile,exist=alive)
    if (.not. alive) then
        write(*,*) "Navigation flie: """//trim(NavFile)//""" doesn't exist!"
        pause
        stop
    end if
    
    NavID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=NavID, file=NavFile,action="read")
    
    do while(.true.)
        read(NavID,fmt="(A)",iostat=error) line
        if (error /=0) exit
        keyword=line(61:80)
        if (index(keyword,"RINEX VERSION / TYPE") /= 0) then
            read(line,fmt="(5X,I1)") NavHead%Version
        elseif (index(keyword,"ION ALPHA") /= 0) then
            read(line,fmt="(2X4D12.4)") NavHead%Alpha
        elseif (index(keyword,"ION BETA") /= 0) then
            read(line,fmt="(2X4D12.4)") NavHead%Beta
        elseif ((index(keyword,"IONOSPHERIC CORR") /= 0) .and. (index(line,"GPSA") /= 0)) then
            read(line,fmt="(5X4D12.4)") NavHead%Alpha   ! This is for RINEX 3.XX
        elseif ((index(keyword,"IONOSPHERIC CORR") /= 0) .and. (index(line,"GPSB") /= 0)) then
            read(line,fmt="(5X4D12.4)") NavHead%Beta
        elseif (index(keyword,"END OF HEADER") /= 0) then
            exit
        end if
    end do
    return
end subroutine