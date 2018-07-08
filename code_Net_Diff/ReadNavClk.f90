

subroutine ReadNavClk(NavFile)
use MOD_FileID
use MOD_NavClk
    implicit none
     character(200) :: NavFile
    logical :: alive
    integer :: error
    character(100) :: line

    integer :: week, prn, prn0, i, temp
    real(8) :: sow, t0, t1

     inquire(file=NavFile,exist=alive)
    if (.not. alive) then
        write(*,*) "COMPASS Navigation clock flie: """//trim(NavFile)//""" doesn't exist!"
        pause
        stop
    end if

    NavID_C=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=NavID_C, file=NavFile, action="read")
    
    i=0
    prn0=0
    do while(.true.)
        read(NavID_C,fmt="(I2)",iostat=error,end=200, advance='no') prn
        if (error /=0) exit
        if (PRN /=PRN0) then
            prn0=prn
            i=0
        end if
        i=i+1
         read(NavID_C,fmt=*,iostat=error,end=200) temp, temp, temp, &
                 & NavClk(PRN)%week(i), NavClk(PRN)%sow(i), NavClk(PRN)%a0(i), NavClk(PRN)%a1(i)
         NavClk(PRN)%week(i)=NavClk(PRN)%week(i)+1356
    end do

    200 return

end subroutine