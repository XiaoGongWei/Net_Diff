!=============  UTC2GPST  ===========
!
! PURPOSE:
!     Transfer from UTC to GPST
!
! REFERENCE:
!     http://www.navipedia.net/index.php/Julian_Date
! 
! INPUTS:
!       UTC                 [yy(yyyy),mon,d,h,min,s]
! 
! OUTPUTS:
!       GPSweek          GPS week
!       GPSsec             GPS second

! Written by: Yize Zhang
!===========  End of Header  ==========

subroutine UTC2GPST(year, m, d,  h, min, s, GPSweek, GPSsec)
implicit none
    ! Intent in
    integer ::  year, m, d, h, min
    real(8) ::  s
    ! Intent out
    integer, intent(out) :: GPSweek
    real(8), intent(out) :: GPSsec
    ! Local variable
    integer(2) :: y, mon
    real(8) :: MJD
    
    y=year
    mon=m
    if (y<80) then
        y=y+2000
    else if ((y>=80) .and. (y<100)) then
        y=y+1900
    end if
    
    if (mon<=2) then
        y=y-1
        mon=mon+12
    end if
    
    MJD=int(365.25d0*y)+int(30.6001d0*dble(mon+1))+d+(dble(h)+(dble(min)+s/60.0d0)/60.0d0)/24.0d0+1720981.5d0-2400000.5d0;
    GPSweek=int((MJD-44244.0d0)/7.0d0);
!    GPSsec=int(mod((MJD-44244.0d0),7.0d0))*86400.0d0+dble(mod(h,24))*3600.0d0+dble(min)*60.0d0+s;
    GPSsec=nint((MJD-44244.0d0-GPSweek*7.d0)*86400.d0*1000.d0)/1000.d0   ! The precision is 0.001s
end subroutine