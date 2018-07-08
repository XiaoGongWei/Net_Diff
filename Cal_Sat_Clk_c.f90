! ================  Cal_Sat_Clk_c  ===============
! PURPOSE:
!      Calculate the BEIDOU satellite clcok using the navigation data
!
! INPUTS:
!    GPSweek           GPS week
!    GPSsec              GPS second
!    PRN                   satellite number, integer
!
! OUTPUT:
!    Sat_Clk              satellite clock, unit in second
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ==================  End of Header  ===============

subroutine Cal_Sat_Clk_c(GPSweek, GPSsec, PRN, Sat_Clk,toe)
use MOD_constant
use MOD_NavData
use MOD_VAR
implicit none
    ! Intent in
    integer :: PRN
    integer :: GPSweek
    real(8) :: GPSsec
    ! Intent out
    real(8) ::  Sat_Clk
    
    type(type_NavPRN) :: tempNav
    real(8) :: t, tkr, tk
    real(8) :: URE(14),toe
!    real(8) :: e, A, n0,n
!    real(8) :: Mk, Ek0, Ek, vk, PHI_k
!    real(8) :: duk, dik, drk
!    real(8) :: uk, radius, ik
!    real(8) :: x, y, z, Omega
!    real(8) :: F, Rotation(3,3), t2
    integer :: i
    real(8) :: delaytmp
    
    if (clktype==1) then
        delaytmp=delay
    elseif (clktype==2) then
        delaytmp=delay-1620.d0  ! omc切换星历时刻，33min
    end if
    tkr=(GPSweek-NavData_C(PRN)%Nav(1)%GPSweek)*604800.d0+GPSsec-NavData_C(PRN)%Nav(1)%GPSsec
    do i=2,Ubound(NavData_C(PRN)%Nav,dim=1)
        t=(GPSweek-NavData_C(PRN)%Nav(i)%GPSweek)*604800.d0+GPSsec-NavData_C(PRN)%Nav(i)%GPSsec
        if ((t-delaytmp>=-0.2d0)  .and. (abs(t-delaytmp)<=abs(tkr)) .and. (NavData_C(PRN)%Nav(i)%Health==0.d0) ) then
            if ((proc_mod/=5) .and. (NavData_C(PRN)%Nav(i)%IODC<10.d0))  cycle
            tkr=t  ! 取最新的星历,t扣除了接收机钟差
            tempNav=NavData_C(PRN)%Nav(i)
        elseif (t<-0.2d0) then
            exit
        end if
!        if (dabs(t)<=dabs(tkr)) then
!            tkr=t ! 取距离最近时刻的星历
!            tempNav=NavData_C(PRN)%Nav(i)
!        else if ( (dabs(t) - dabs(tkr))>700.d0 ) then
!            exit
!        end if
    end do
    toe=tempNav%GPSsec
    if ( (dabs(tkr-delaytmp)>3600.2d0) ) then
        Sat_Clk=9999.d0
        return
    end if

     ! Satellite clock
    tk=tkr
    Sat_Clk=tempNav%a0+tk*tempNav%a1+tk*tk*tempNav%a2

    ! 北斗广播星历卫星钟TGD改正
    if (index(ObsCombine,"PC")/=0) then
        if (freq_comb=='L1L2') then
            Sat_Clk=Sat_Clk - 2.4872*tempNav%TGD(1) + 1.4872*tempNav%TGD(2)   ! B3==>B1B2  f1^2/(f1^2-f2^2)*T13+f2^2/(f1^2-f2^2)*[T23]
            Sat_Clk=Sat_Clk + ( 2.4872*BEB(1,PRN)- 1.4872*BEB(2,PRN) )/c    ! Broadcast Ephemeris Bias
        elseif (freq_comb=='L1L3') then
            Sat_Clk=Sat_Clk - 2.94368*tempNav%TGD(1)                    ! B3==>B1B3 f1^2/(f1^2-f3^2)*T13
            Sat_Clk=Sat_Clk + ( 2.94368*BEB(1,PRN)- 1.94368*BEB(3,PRN) )/c    ! Broadcast Ephemeris Bias
        else  if (freq_comb=='L2L3') then
            Sat_Clk=Sat_Clk + 9.589*tempNav%TGD(2)        ! B3==>B2B3 f2^2/(f2^2-f3^2)*T23
            Sat_Clk=Sat_Clk + (- 9.589*BEB(2,PRN) + 10.5895*BEB(3,PRN) )/c    ! Broadcast Ephemeris Bias
        end if
    elseif ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0))then   
        Sat_Clk=Sat_Clk - tempNav%TGD(1)                       ! B3==>B1
        Sat_Clk=Sat_Clk + BEB(1,PRN)/c    ! Broadcast Ephemeris Bias
    elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) then
        Sat_Clk=Sat_Clk - tempNav%TGD(2)                      ! B3==>B2
        Sat_Clk=Sat_Clk + BEB(2,PRN)/c    ! Broadcast Ephemeris Bias
    elseif ((index(ObsCombine,"P3")/=0) .or. (index(ObsCombine,"G3")/=0)) then
        Sat_Clk=Sat_Clk + BEB(3,PRN)/c    ! Broadcast Ephemeris Bias
    end if

!     !  X31a卫星钟TGD改正
!    if (freq_comb=='L1L2') then
!        Sat_Clk=Sat_Clk - 1.4872*TGD(1, PRN) - TGD(2,PRN)   ! B3==>B1B2  f1^2/(f1^2-f2^2)*T13+f2^2/(f1^2-f2^2)*[T13-T12]
!    elseif (freq_comb=='L1L3') then
!        Sat_Clk=Sat_Clk - 2.94368*TGD(2, PRN)                    ! B3==>B1B3 f1^2/(f1^2-f3^2)*T13
!    else  if (freq_comb=='L2L3') then
!        Sat_Clk=Sat_Clk + 9.589*(TGD(2, PRN) - TGD(1,PRN))  ! B3==>B2B3 f2^2/(f1^2-f2^2)*T12+f2^2/(f2^2-f3^2)*T12-f3^2/(f2^2-f3^2)*T13
!    end if

end subroutine