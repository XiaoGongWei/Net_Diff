! ============== Get_BEB =================
!
! PURPOSE:
!          Get BEBs of BeiDou Bradcast Ephemeris Bias BEB file.
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================  End of header  =================

subroutine Get_BEB
use MOD_VAR
use MOD_FileID
use MOD_FileDir
implicit none
    logical   :: alive
    character(200) :: line
    integer   :: error,PRN
    real(8)   :: temp(3)

    
    BEB=0.d0
    inquire(file=BEBFile,exist=alive)
    if (.not. alive) then
!        write(*,*) "WARN: BEB file doesn't exist, BEB of BeiDou will default as 0!"
        return
    end if
    BEBID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=BEBID, file=BEBFile)
    
    do while(.true.)
        read(unit=BEBID,fmt='(A)', iostat=error) line
        if (error /=0) exit
        read(line,*) PRN,temp(1:3)
        if (PRN>35) cycle
        BEB(1,PRN)=temp(1)
        BEB(2,PRN)=temp(2)
        BEB(3,PRN)=temp(3)
    end do

    close(BEBID)

end subroutine

