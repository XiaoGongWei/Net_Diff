! ================  Cal_Sat_Clk_n  ===============
! PURPOSE:
!      Calculate the GPS satellite clcok using the navigation data
!
! INPUTS:
!    System               Sytem flag, character(1)
!    GPSweek           GPS week
!    GPSsec              GPS second
!    PRN                   satellite number, integer
!
! OUTPUT:
!    Sat_Clk              satellite clock, unit in second
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ==================  End of Header  ===============

subroutine Cal_Sat_Clk_n(System, GPSweek, GPSsec, PRN, Sat_Clk)
use MOD_constant
use MOD_NavData
use MOD_VAR
implicit none
    ! Intent in
    character(1) :: System
    integer PRN
    integer :: GPSweek
    real(8) :: GPSsec
    ! Intent out
    real(8) ::  Sat_Clk
    
    type(type_NavPRN) :: tempNav
    real(8) :: t, tkr, tk
!    real(8) :: e, A, n0,n
!    real(8) :: Mk, Ek0, Ek, vk, PHI_k
!    real(8) :: duk, dik, drk
!    real(8) :: uk, radius, ik
!    real(8) :: x, y, Omega
!    real(8) :: F, Rotation(3,3), t2
    integer :: i

    if (System=='G') then
        tkr=(GPSweek-Navdata(PRN)%Nav(1)%GPSweek)*604800.d0+GPSsec-NavData(PRN)%Nav(1)%GPSsec
        do i=2,Ubound(NavData(PRN)%Nav,dim=1)
            t=(GPSweek-Navdata(PRN)%Nav(i)%GPSweek)*604800.d0+GPSsec-NavData(PRN)%Nav(i)%GPSsec
            if (dabs(t)<=dabs(tkr)-0.1d0 .and. ((NavData(PRN)%Nav(i)%Health==0.d0) .or. (IsCNAV))) then
                tkr=t    ! 取距离最近时刻的星历
                tempNav=NavData(PRN)%Nav(i)
            else if ( (dabs(t) - dabs(tkr))>700.d0 ) then
                exit
            end if
!            if ((t>=-0.1d0)  .and. (abs(t)-abs(tkr)<20.d0) .and. ((NavData(PRN)%Nav(i)%Health==0.d0) .or. (IsCNAV)) ) then
!                tkr=t   ! 取最新的星历,t扣除了接收机钟差，排除59‘44’‘的星历
!                tempNav=NavData(PRN)%Nav(i)
!            elseif (t<-0.2d0) then
!                exit
!            end if
        end do
        if ((dabs(tkr)>7220.d0) ) then
            Sat_Clk=9999.d0
            return
        end if
    elseif (System=='J') then  ! QZSS
        tkr=(GPSweek-NavData_J(PRN)%Nav(1)%GPSweek)*604800.d0+GPSsec-NavData_J(PRN)%Nav(1)%GPSsec
        do i=2,Ubound(NavData_J(PRN)%Nav,dim=1)
            t=(GPSweek-NavData_J(PRN)%Nav(i)%GPSweek)*604800.d0+GPSsec-NavData_J(PRN)%Nav(i)%GPSsec
    !        if (dabs(t)<=dabs(tkr)-0.1d0) then
    !            tkr=t    ! 取距离最近时刻的星历
    !            tempNav=NavData_J(PRN)%Nav(i)
    !        else if ( (dabs(t) - dabs(tkr))>700.d0 ) then
    !            exit
    !        end if
            if ((t>=-0.2d0)  .and. (abs(t)-abs(tkr)<20.d0) ) then !  .and. (NavData_J(PRN)%Nav(i)%Health==0.d0)
                tkr=t   ! 取最新的星历,t扣除了接收机钟差，排除59‘44’‘的星历
                tempNav=NavData_J(PRN)%Nav(i)
            elseif (t<-0.2d0) then
                exit
            end if
        end do
        if ((dabs(tkr)>3620.d0) ) then
            Sat_Clk=9999.d0
            return
        end if
    elseif (System=='I') then  ! IRNSS
        tkr=(GPSweek-NavData_I(PRN)%Nav(1)%GPSweek)*604800.d0+GPSsec-NavData_I(PRN)%Nav(1)%GPSsec
        do i=2,Ubound(NavData_I(PRN)%Nav,dim=1)
            t=(GPSweek-NavData_I(PRN)%Nav(i)%GPSweek)*604800.d0+GPSsec-NavData_I(PRN)%Nav(i)%GPSsec
    !        if (dabs(t)<=dabs(tkr)-0.1d0) then
    !            tkr=t    ! 取距离最近时刻的星历
    !            tempNav=NavData_I(PRN)%Nav(i)
    !        else if ( (dabs(t) - dabs(tkr))>700.d0 ) then
    !            exit
    !        end if
            if ((t>=-0.2d0)  .and. (abs(t)-abs(tkr)<20.d0) .and. (NavData_I(PRN)%Nav(i)%Health==0.d0) ) then 
                tkr=t
                tempNav=NavData_I(PRN)%Nav(i)
            elseif (t<-0.2d0) then
                exit
            end if
        end do
        if ((dabs(tkr)>7220.d0) ) then
            Sat_Clk=9999.d0
            return
        end if
    end if

    ! Satellite clock
    tk=tkr
    Sat_Clk=tempNav%a0+tk*tempNav%a1+tk*tk*tempNav%a2

    ! 广播星历卫星钟TGD改正 ！ Reference: GPS ICD
    ! For IRNSS, L5 is treated as P1 here. But according to IRNSS ICD, S band is P1.
    ! However, SPP residuals and results show that L5 as P1 is better. Need further investigation.
    if ((System/='G') .or. (.not.(IsCNAV))) then ! TGD of CNAV is already corrected in observations
        if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0))then
            Sat_Clk=Sat_Clk - tempNav%TGD(1)                       ! B1B2==>B1
        elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) then
            Sat_Clk=Sat_Clk - (f1/f2)**2*tempNav%TGD(1)                      ! B1B2==>B2
        end if
    end if

end subroutine