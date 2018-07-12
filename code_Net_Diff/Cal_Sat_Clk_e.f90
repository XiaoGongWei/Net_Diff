! ================  Cal_Sat_Clk_e  ===============
! PURPOSE:
!      Calculate the GALILEO satellite clcok using the navigation data
!
! REFERENCE:
!    http://www.navipedia.net/index.php/Galileo_Navigation_Message
!    Galileo SIS ICD, 2010
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

subroutine Cal_Sat_Clk_e(GPSweek, GPSsec, PRN, Sat_Clk)
use MOD_constant
use MOD_NavData
use MOD_VAR
implicit none
    ! Intent in
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

    tkr=(GPSweek-Navdata_E(PRN)%Nav(1)%GPSweek)*604800.d0+GPSsec-Navdata_E(PRN)%Nav(1)%GPSsec
    do i=2,Ubound(Navdata_E(PRN)%Nav,dim=1)
        t=(GPSweek-Navdata_E(PRN)%Nav(i)%GPSweek)*604800.d0+GPSsec-Navdata_E(PRN)%Nav(i)%GPSsec
        if (dabs(t)<=dabs(tkr)-0.1d0 .and. (NavData_E(PRN)%Nav(i)%Health==0.d0)) then
            tkr=t    ! 取距离最近时刻的星历
            tempNav=Navdata_E(PRN)%Nav(i)
        else if ( (dabs(t) - dabs(tkr))>700.d0 ) then
            exit
        end if
!        if ((t>=-0.1d0)  .and. (abs(t)-abs(tkr)<20.d0) .and. (NavData_E(PRN)%Nav(i)%Health==0.d0) ) then
!            tkr=t   ! 取最新的星历,t扣除了接收机钟差，! .and. (mod(NavData_E(PRN)%Nav(i)%Code,2.d0)==0.d0) Code=258 means F/NAV, E1/E5a; Code =513 or 517 mens I/Nav, E1/E5b
!            tempNav=Navdata_E(PRN)%Nav(i)
!            if ((dabs(tkr)<3620.d0)  .and. (tempNav%Health==0.d0)  .and. (tempNav%a0/=0.d0) ) then
!                exit
!            end if
!!        elseif (t<-0.2d0) then
!!            exit
!        end if
    end do
    if ((dabs(tkr)>3620.d0) ) then
        Sat_Clk=9999.d0
        return
    end if

    ! Satellite clock
    tk=tkr
    Sat_Clk=tempNav%a0+tk*tempNav%a1+tk*tk*tempNav%a2
    
    ! BGD correction ！ Reference: Galileo ICD 5.1.3~5.1.5
    ! F/Nav refered to E1/E5a BGD  
    ! I/Nav refered to E1/E5b BGD  E1/E5a  E1/E5b
    ! In fact, differences between the E1–E5a and E1–E5b biases are
    ! sufficiently small (maybe less thant 0.3ns) 
    !      Montenbruck O, Hauschild A, Steigenberger P (2014) Differential
    ! code bias estimation using multi-GNSS observations and global
    ! ionosphere maps. In: Proceedings of ION ITM 2014, San Diego,CA)
    ! Notice: IGS Precise clock refered to E1E5a
    ! The following clock is based on F/Nav, i.e. refered to E1/E5a.
    if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
!        if (tempNav%TGD(2)==0.d0)  tempNav%TGD(2)=tempNav%TGD(1)  ! If no BGD(E1,E5b), BGD(E1,E5a) insteaded
        Sat_Clk=Sat_Clk - tempNav%TGD(2)                       ! E1/E5a==>E1
    elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) then
        Sat_Clk=Sat_Clk - (f1/f2)**2*tempNav%TGD(1)     ! E1/E5a==>E5a, assume  that clock in I/Nav=F/nav
    elseif ((index(ObsCombine,"P3")/=0) .or. (index(ObsCombine,"G3")/=0)) then
!        if (tempNav%TGD(2)==0.d0)  tempNav%TGD(2)=tempNav%TGD(1)  ! If no BGD(E1,E5b), BGD(E1,E5a) insteaded
        Sat_Clk=Sat_Clk - (f1/f2)**2*tempNav%TGD(2)     ! E1/E5a==>E5b
    end if

end subroutine