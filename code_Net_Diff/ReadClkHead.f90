! =======ReadClkHead=====

! Read the Header of clk File

! Written by Yize Zhang
! ======End of header=======

subroutine ReadClkHead(ClkFile)
use MOD_FileID
implicit none
    character(200) ClkFile
    logical :: alive
    integer :: error
    character(100) :: line
    
    inquire(file=ClkFile,exist=alive)
    if (.not. alive) then
        write(*,*) "Clock file: """//trim(ClkFile)//""" doesn't exist!"
        pause
        stop
    end if
    
    ClkID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=ClkID, file=ClkFile)
    
    do while(.true.)
        read(unit=ClkID,fmt="(A)") line
        if (index(line,"END OF HEADER") /= 0) then
            exit
        end if
    end do
    return
end subroutine