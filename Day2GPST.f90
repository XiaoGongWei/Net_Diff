!==========   From Day to GPST  ==========
!
! PURPOSE:
!     Change year and doy to GPS time
! 
! NOTE:
!      This function is useful during year 1901-2099
!
! INPUTS:
!       year
!       doy
!
! OUTPUTS:
!       GPSweek              GPS week
!       GPSday                 GPS day
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!=========== End of Header ==========

subroutine Day2GPST(year, doy, GPSweek, GPSday)
implicit none
    ! Intent in
    integer ::  year, doy

    ! Intent out
    integer, intent(out) :: GPSweek, GPSday
    ! Local variable
    integer(2) :: y
    integer :: nday
    
    y=year
    if (y<80) then
        y=y+2000
    else if ((y>=80) .and. (y<100)) then
        y=y+1900
    end if
    
!    ***Tot. nr of days since 05.01.1980*******************************
      nday=(year-1980)*365+doy-5
 !     **Correction for leap-years***************************************
      nday=nday+(year-1980)/4
      IF (MOD((year-1980),4)==0) nday=nday-1 
      GPSweek = INT(nday/7)
      GPSday=nday-GPSweek*7
!    MJD=int(365.25d0*y)+doy-679019
!    GPSweek=int((MJD-44244)/7)
!    GPSday=int(mod((MJD-44244.0d0),7.0d0))
    
end subroutine