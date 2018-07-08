! ============== Get_GLO_Fre ============
!
! PURPOSE:
!     This subroutine is used to read the GLONASS 
! frequency information and leap seconds.
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!   ============ End of Header ==========

subroutine Get_GLO_Fre
use MOD_Var
use MOD_FileID
use MOD_FileDir
use MOD_GLO_Fre
implicit none
    ! Local cariables
    logical alive
    character(300) :: line
    integer :: year, mon, day, doy
    real(8) :: MJD0, MJD, seconds
!    integer :: PRN, day_st, day_end, Channel, 
    
    call UTC2MJD(int_year, 1, int_doy, 0, 0, 0.d0, MJD)
    inquire(file=GLOFreFile,exist=alive)
    if (.not. alive) then
        write(*,*) "GLONASS frequency file: """//trim(GLoFreFile)//""" doesn't exist!"
        pause
        stop
    end if
    GLOFreID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=GLOFreID, file=GLOFreFile)
    
    do while(.true.)
        read(unit=GLOFreID,fmt='(A)',end=200) line
        if (index(line,"++GLONASS Information") /= 0)  then
            do while(.true.)
                read(unit=GLOFreID,fmt='(A)',end=200) line
                if (len_trim(line)==0) cycle
                if (line(1:1)=="*") cycle
                if (line(1:1)=="-") exit
                read(line,*) year, mon, day
                READ(GLOFreID,*) Fre_Chann
                call UTC2MJD(year,mon,day, 0, 0, 0.d0, MJD0)
                if (MJD0>MJD) exit
!                read(line(3:4),*) PRN
!                read(line(18:24),*) day_st
!                read(line(32:38),*) day_end
!                read(line(55:56),*) Channel
!                if (day_end==0) day_end=9999999
!                if ( (day>=day_st) .and. (day<=day_end) ) then
!                    Fre_Chann(PRN)=Channel
!                end if
            end do
        elseif (index(line,"++Leap seconds") /= 0)  then
            do while(.true.)
                read(unit=GLOFreID,fmt='(A)',end=200) line
                if (line(1:1)=="*") cycle
                if (line(1:1)=="-") exit
                read(line,*) year, doy, seconds
                if ( (year>=int_year) .and. (doy>int_doy) ) exit
                Leap_sec=seconds
            end do
        end if
    end do
    200 close(GLOFreID)
    if (Leap_sec==0.d0) then
        write(*,*) "Leap_sec not found in day ", day
        pause
        stop
    end if
    return
    
end subroutine
