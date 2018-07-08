! ==============NavCorr =============
! Purpose:
!           get navgation corrctions
!
! Inputs:
!      PRN
!      sow               seconds of week
!      LoP               the length of predicte ephemeris
!      Sat_Coor(3)
!      Sat_Clk
!
! Output:
!      Sat_Coor(3)
!       Sat_Clk
!
! Written by Yize Zhang
! ===============================

subroutine NavCorr(PRN, sow, LOP, Sat_Coor, Sat_Clk)
use MOD_FileID
implicit none
    ! Intents:
    integer  :: PRN
    real(8)   :: sow
    real(8)   :: LOP
    real(8)   :: Sat_Coor(3)
    real(8)   :: Sat_Clk
    ! Local variables
    integer :: i, error
    character(50) :: line
    logical :: flag
    real(8) :: sod, t
    real(8) :: sec, coe_X(3), coe_Y(3), coe_Z(3), coe_Clk(3)
    real(8) :: dXYZ(3), dClk

    flag=.false.
    sod=mod(sow,86400.d0)  ! second of day
    do while (.true.)
        read(modelID(PRN),*, iostat=error) line, sec
        if (error /=0) exit
        t=sod-sec
        if (t>LOP+0.1d0) then   ! 超过模型预报时间
            read(modelID(PRN),*) 
            read(modelID(PRN),*)
            read(modelID(PRN),*)
            read(modelID(PRN),*)
        elseif (t+0.1d0>0.d0) then   ! 在预报时间内
            flag=.true.
            read(modelID(PRN),*) line, coe_X(1:3)
            read(modelID(PRN),*) line, coe_Y(1:2)
            read(modelID(PRN),*) line, coe_Z(1:3)
            read(modelID(PRN),*) line, coe_Clk(1:3)
            backspace(modelID(PRN))
            backspace(modelID(PRN))
            backspace(modelID(PRN))
            backspace(modelID(PRN))
            backspace(modelID(PRN))
            exit
        else                   ! 还没到预报时间
            backspace(modelID(PRN))
            exit
        end if
    end do

    if (flag) then
        t=t/3600.d0
        dXYZ(1)=DOT_PRODUCT(coe_X, (/1.d0, t, t*t /) )
        dXYZ(2)=DOT_PRODUCT(coe_Y, (/1.d0, t, t*t /) )
        dXYZ(3)=DOT_PRODUCT(coe_Z, (/1.d0, t, t*t /) )
        dClk=DOT_PRODUCT(coe_Clk, (/1.d0, t, t*t /) )
        Sat_Coor=Sat_Coor - dXYZ
        Sat_Clk=Sat_Clk - dClk/299792458.0d0
    else
        Sat_Coor=9999.d0
        Sat_Clk=9999.d0
    end if

    return
end subroutine