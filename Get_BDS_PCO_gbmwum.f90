! ========== Get_BDS_PCO_gbmwum ========
!
! PURPOSE:
!      GBM and WUM MGEX products use different PCO for BeiDou.

! Reference:
!       ftp://ftp.gfz-potsdam.de/GNSS/products/mgex/SISRE/
!      郭靖. 姿态、光压和函数模型对导航卫星精密定轨影响的研究[D]. 武汉大学, 2014.

! INPUTS:
!      year             year
!      doy              doy
!      ac                analysis center type, 1 stands for gbm and 2 stands for wum
!
! OUTPUT(as a public module data):
!      Ant               A type data, See MOD.F90
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!   ============ End of Header ==========

subroutine Get_BDS_PCO_gbmwum(year, doy, ac)
use MOD_Ant
    integer :: year, doy, ac
    integer :: PRN
    real(8) :: PCO(14,3,2)=0.d0

    PCO(6:14,1,1)=0.549d0
    PCO(6:14,3,1)=(/3.049, 3.237, 3.843, 3.974, 3.882, 2.069, 2.313, 3.882, 2.312/)
    PCO(6:10,1,2)=0.586d0
    PCO(11:14,1,2)=0.575d0
    PCO(6:14,3,2)=(/2.514, 2.722,  3.440, 3.552, 4.087, 1.991, 2.249, 2.896, 2.144/)
    if (year*100+doy<2016285) then
        PCO(13,3,1)=2.202
        PCO(13,3,2)=2.026
    elseif (year*100+doy<2017253) then
        PCO(13,1:3,2)=(/0.6, 0.0, 1.1/)
    end if

    do PRN=6,14
        Ant(PRN)%PCO(:,1)=PCO(PRN, :, ac)
        Ant(PRN)%PCO(:,2)=PCO(PRN, :, ac)
    end do


end subroutine