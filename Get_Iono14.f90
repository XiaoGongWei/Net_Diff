! ================================
! PURPOSE:
!     Get ionospheric 14 parameter.
!
! ========End of Header============

subroutine Get_Iono14(Iono14file)
use MOD_VAR
use MOD_FileID
implicit none
    character(200) :: Iono14file, line
    integer   :: error
    logical :: alive
    integer(1) :: t
    real(8) :: temp(1:14)

    inquire(file=Iono14file,exist=alive)
    if (.not. alive) then
        write(*,*) "Iono_14_para file: """//trim(Iono14file)//""" doesn't exist!"
        pause
        stop
    end if
    IonID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=IonID, file=Iono14file)
    
    do while(.true.)
        read(unit=IonID,fmt='(A)',iostat=error) line
        if (error /=0) exit
        read(line,*) t,temp(1:14)
        t=int(t/2)+1
        iono_14para_day(t,1:2)=temp(13:14)
        iono_14para_day(t,3:14)=temp(1:12)
    end do

    return
end subroutine