! =================  Ocean Load  =================
! Purpose: 
!      Calculate the ocean load.
!
! REFERENCE:
!      http://www.navipedia.net/index.php/Ocean_loading
! 
! Adapted from: A_RTK, tide_disp_oloas.f90(by Junping Chen),
!     which is referenced from EPOS.
!
! Inputs:
!     Flag          CMC flag, if Flag='CMC', calculate the CMC.
!     IYEAR       Year(2 digit)
!     DAY          Day of year, including the fraction of the day
!     OLC          Ocean Load Coefficient
!
! Outputs:
!    dNEU        Change of ocean load in North-East-Up
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================  End of Header  ==============

subroutine OceanLoad(Flag,IYear,Day,OLC, dNEU)
implicit none
    ! Intent in
    character(3) :: Flag
    integer :: IYEAR
    real(8) :: Day
    real(8) :: OLC(11,6)
    ! Intent out
    real(8) :: dNEU(3)
    ! Local variables
    real(8) :: Angle(11)
    integer :: i
    
    dNEU(1:3) = 0.0d0
    ! for station ENU components, for CMC xyz in earth-fixed system
    if (Flag=="CMC" .or. Flag=='cmc') then
        !call 
    end if
    
    ! ocearg expects UT1, but use UTC instead (max error 1 sec = 10 ppm)
    call ARG2(IYear,Day,Angle)
    do i=1,11
       dNEU(2)=dNEU(2)+OLC(i,1)*(dcos(angle(i)-OLC(i,4)))   ! In A_RTK, it is E-N-U,we change it to NEU
       dNEU(1)=dNEU(1)+OLC(i,2)*(dcos(angle(i)-OLC(i,5)))
       dNEU(3)=dNEU(3)+OLC(i,3)*(dcos(angle(i)-OLC(i,6)))
    end do
    return
end subroutine