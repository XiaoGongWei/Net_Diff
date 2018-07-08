! =======Get EOP(Earth Orientation Parameters)========
!
! PURPOSE:
!     Read the file 'EOP.txt' and get the Earth Orientation Parameters.
!
! INPUT:
!    MJD       usually in the middle of the day
!
! OUTPUTS:
!    RMJD
!    X           unit the same as in EOP.txt
!    Y
!    dUT1    UT1-UTC
!    dX
!    dY
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!   =============== End of Header ==============

subroutine Get_EOP(MJD, RMJD,X, Y, dUT1, dX, dY)
use MOD_FileID
use MOD_FileDir
implicit none
    ! Intent in
    real(8) :: MJD
    ! Intent out:
    real(8) :: RMJD(2), X(2), Y(2), dUT1(2), dX(2), dY(2)
    
    ! Local variables
!    character(10) :: EOPFile='EOP.txt'
    logical :: alive
    character(100) :: line
    real(8) :: tempMJD
    
    inquire(file=EOPFile,exist=alive)
    if (.not. alive) then
        write(*,*) "EOP file: """//trim(EOPFile)//""" doesn't exist, EOP will be set as 0."
        pause
    end if
    EOPID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=EOPID, file=EOPFile)
    
    do while(.true.)
        read(unit=EOPID,fmt='(A)',end=200) line
        if ( line(1:5)=='-----') exit
    end do
    
    do while (.true.)
        read(unit=EOPID,fmt='(12X,F7.0, 2F11.6, F12.7, 12X, 2F11.6)',end =300) RMJD(2), X(2), Y(2), dUT1(2), dX(2), dY(2)
        if ((MJD-RMJD(2)<=1.d0) .and. (MJD-RMJD(2)>=0.d0) )then
            RMJD(1)=RMJD(2)
            X(1)=X(2)
            Y(1)=Y(2)
            dUT1(1)=dUT1(2)
            dX(1)=dX(2)
            dY(1)=dY(2)
            read(unit=EOPID,fmt='(12X,F7.0, 2F11.6, F12.7,12X,2F11.6)',end=300) RMJD(2), X(2), Y(2), dUT1(2), dX(2), dY(2)
            goto 300
        else if (MJD<RMJD(2)) then
!            write(*,*) "Current MJD is earlier than the first MJD in EOP file, EOP will be set as 0."
!            pause
        end if
    end do
    
    
    200 close(EOPID)
    write(*,*) "EOP file format error! Start of ""-------- ""not found! "
    pause
    stop
    300 close(EOPID)
    if (MJD>RMJD(2)) then
        write(*,*) "Current MJD is later than the last MJD in EOP file, EOP will be set as 0."
        !pause
    end if
    return
    
end subroutine