! =============== PCOR =========================
! PURPOSE:
!          Get PCOR
!
! ==============================================
module MOD_ClkSmooth
implicit none
   type type_PCOR
        integer :: epoch=0
        integer :: st=1
        real(8) :: sow(2880)=0.d0
        real(8) :: PCOR(2880)=0.d0
   end type
    type type_ClkSmooth
        logical :: flag=.false.
        type(type_PCOR) :: val(32)
    end type
    type(type_ClkSmooth) :: ClkSmooth
end module

subroutine PCOR(GPSweek, GPSsec, PRN, Sat_Clk)
use MOD_constant
use MOD_VAR
use MOD_ClkSmooth
implicit none
    ! Intent in
    integer PRN
    integer :: GPSweek
    real(8) :: GPSsec
    ! Intent out
    real(8) ::  Sat_Clk
    ! Local variables
    integer :: week, i, satid
    real(8) :: sow, t,temp
    logical :: alive

    if (ClkSmooth%flag==.false.) then
        open(unit=22, file='J:\Â½Ì¬ÍøÊý¾Ý\clksmooth',action='read')
        do while (.true.)
            read(unit=22,fmt='(I4,I7,F7.0,F12.3,2F13.3)',end=100) SatID, week, sow, temp, temp, t
            ClkSmooth%val(SatID)%epoch=ClkSmooth%val(SatID)%epoch+1
            ClkSmooth%val(SatID)%sow(ClkSmooth%val(SatID)%epoch)=sow
            ClkSmooth%val(SatID)%PCOR(ClkSmooth%val(SatID)%epoch)=t
        end do
        100 ClkSmooth%flag=.true.
        close(unit=22)
    end if
    
    do i=ClkSmooth%val(PRN)%st,ClkSmooth%val(PRN)%epoch
        if (dabs(GPSsec-ClkSmooth%val(PRN)%sow(i))<1.d0) then
            Sat_Clk=Sat_Clk+ClkSmooth%val(PRN)%PCOR(i)/c
            ClkSmooth%val(PRN)%st=i
            exit
        elseif (ClkSmooth%val(PRN)%sow(i)-GPSsec>Interval) then
            exit
        end if
    end do

    return
end subroutine