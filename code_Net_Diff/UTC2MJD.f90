!======From UTC to MJD======

! Input:
!       UTC                 [yy(yyyy),mon,d,h,min,s]
! Outputs:
!       MJD                

! Written by: Yize Zhang
!===========End of Header==========

subroutine UTC2MJD(year,m,d, h, min, s, MJD)
implicit none
    ! Intent in
    integer ::  year, m, d, h, min
    real(8) ::  s
    ! Intent out
    real(8) :: MJD
    ! Local variable
    integer(2) :: y,mon

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
    
    MJD=int(365.25d0*y)+int(30.6001d0*dble(mon+1))+d+(dble(h)+(dble(min)+s/60.0d0)/60.0d0)/24.0d0+1720981.5d0-2400000.5d0
    
end subroutine