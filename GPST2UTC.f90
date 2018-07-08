!============ From GPST to UTC ===============
!
! PURPOSE:
!        Transfer GPST to UTC
!
! INPUTS:
!       GPST               week, sow in GPST           
! 
! OUTPUTS:
!       MJD
!       UTC               [yyyy,mon,d,h,min,s]           
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!==============  End of Header =============

subroutine GPST2UTC(week, sow, MJD,year,mon,day,hour,min,sec)
implicit none
    ! Intent in:
    integer :: week
    real(8) :: sow
    ! Intent out
    real(8) :: MJD
    integer :: year, mon, day, hour, min, doy
    real(8) :: sec
    ! Local varibles
    integer :: a,b,c,d,e
    REAL(8) :: FD
    
    MJD=week*7.d0+44244.0d0+sow/86400.d0   ! In GPST
    
!    MJD=MJD+1356*7.d0  !+14.d0/86400.d0      !In BDT
    
    a=int(MJD+1.d0+1.d-10)+2400000
    b=a+1537
    c=int(dble(b-122.1d0)/365.25d0)
    d=int(365.25d0*dble(c))
    e=int(dble(b-d)/30.6001d0)
    day=int(b-d-int(30.6001d0*dble(e)))
    mon=e-1-12*int(dble(e)/14.d0+1.d-10)
    year=c-4715-int(dble(7+mon)/10.d0+1.d-10)
    
    FD=MJD-int(MJD)
!    hour=int(FD*24.d0)
!    min=int((FD*24.d0-hour)*60.d0)
    hour=int(mod(sow,86400.d0)/3600.d0)
    min=int(mod(sow,3600.d0)/60.d0)
    sec=mod(sow,86400.d0)-dble(hour*3600+min*60)

    return
end subroutine