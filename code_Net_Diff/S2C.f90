! ====================  S2C  ===================
!
! PURPOSE:
!       Transformation matrix from Satellite-fixed system to 
!   Celestrial reference system
!
!  REFERENCE:
!      http://www.navipedia.net/index.php/Satellite_Antenna_Phase_Centre
! 
! INPUTS:
!        Sun_Coor           Coordinate of the sun, in CRS
!        Sat_Coor            Coordinate of the satellite mass center, in TRS
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ====================  End of header  =====================

subroutine S2C(Sun_Coor, Sat_Coor, Sat_Vel, System, PRN)
use MOD_Ant
use MOD_Rotation
use MOD_constant
use MOD_FileID
use MOD_FileDir
implicit none
    ! Intent
    real(8) :: Sun_Coor(3), Sat_Coor(3)
    ! Local variables
    real(8) :: Sat_Coor_C(3)
    real(8) :: Sun_Sat(3)
    real(8) :: Sat_Vel(3),Sat_Vel_C(3), NormalVec(3), beta, mu, yaw, mu_rate, yaw_rate
    integer :: PRN
    character(1) :: System

            ! Satellite mass center coordinate in CRS
            Sat_Coor_C=MATMUL(Rota_T2C, Sat_Coor) 
             ! ez: Satellite fixed reference sysytem unit vector of z in CRS
             !   unit vector pointing from the Satellite Mass Centre (MC) to the earth's centre
            Rota_S2C(1:3,3)= - Sat_Coor_C/dsqrt(DOT_PRODUCT(Sat_Coor_C,Sat_Coor_C) ) 
            ! ey: Satellite fixed reference sysytem unit vector of x in CRS
            !   the cross product of ez with the unit vector from the satellite to sun
            Sun_Sat=Sun_Coor(1:3)-Sat_Coor_C   ! Sun - Satellite in CRS
            call cross( -Sat_Coor_C, Sun_Sat, Rota_S2C(1:3,2))

            Rota_S2C(1:3,2)= Rota_S2C(1:3,2)/dsqrt(DOT_PRODUCT(Rota_S2C(1:3,2), Rota_S2C(1:3,2)))

            ! ex: Satellite fixed reference sysytem unit vector of y in CRS
            !    ex completes the right-handed system
            call cross(Rota_S2C(1:3,2), Rota_S2C(1:3,3), Rota_S2C(1:3,1))
            
            ! For BeiDou IGSO and MEO
            Sat_Vel_C=MATMUL(Rota_T2C, Sat_Vel)
            call cross(Sat_Coor_C, Sat_Vel_C, NormalVec)
            beta=asind(DOT_PRODUCT(Sun_Coor, NormalVec)/dsqrt(DOT_PRODUCT(NormalVec, NormalVec))/dsqrt(DOT_PRODUCT(Sun_Coor, Sun_Coor)))
            mu=acosd(DOT_PRODUCT(-Sun_Coor*cosd(beta), Sat_Coor_C)/dsqrt(DOT_PRODUCT(Sat_Coor_C, Sat_Coor_C))/dsqrt(DOT_PRODUCT(Sun_Coor*cosd(beta), Sun_Coor*cosd(beta))))
            yaw=atan2d(-tand(beta),sind(mu))
            ! mu_rate=acosd(DOT_PRODUCT(-Sun_Coor*cosd(beta), Sat_Coor_C+Sat_Vel_C) &
            ! /dsqrt(DOT_PRODUCT(Sat_Coor_C+Sat_Vel_C, Sat_Coor_C+Sat_Vel_C))/dsqrt(DOT_PRODUCT(Sun_Coor*cosd(beta), Sun_Coor*cosd(beta))))-mu
            ! yaw_rate=mu_rate*tand(beta)*cosd(mu)/(sind(mu)**2+tand(beta)**2)
!            if (System=='G') then
!                write(CSID,'(5X,A5,A2,I2.2,3F10.3)') 'ATT','G',PRN, beta, mu, yaw
!            elseif (System=='R') then
!                write(CSID,'(5X,A5,A2,I2.2,3F10.3)') 'ATT','R',PRN-GNum, beta, mu, yaw
!            elseif (System=='C') then
!                write(CSID,'(5X,A5,A2,I2.2,3F10.3)') 'ATT','C',PRN-GNum-RNum, beta, mu, yaw
!            end if
            ! for BeiDou GEO, always yaw-fixed, orbit-normal mode
            if ( (System=='C') .and. ((PRN-GNum-RNum<6) .or. (PRN-GNum-RNum==17) .or. (PRN-GNum-RNum==18)) ) then
                call cross( -Sat_Coor_C, Sat_Vel_C, Rota_S2C(1:3,2))   ! Y正交于卫星位置-速度平面
                Rota_S2C(1:3,2)= Rota_S2C(1:3,2)/dsqrt(DOT_PRODUCT(Rota_S2C(1:3,2), Rota_S2C(1:3,2)))
                call cross(Rota_S2C(1:3,2), Rota_S2C(1:3,3), Rota_S2C(1:3,1))
            elseif ( (System=='C') .and. (abs(beta)<4.d0) .and.  ((index(SP3File(len_trim(SP3File)-11:len_trim(SP3File)),'gbm')/=0) .or.  &
            (index(SP3File(len_trim(SP3File)-11:len_trim(SP3File)),'wum')==0)) ) then ! for BeiDou IGSO/MEO, yaw-fixed, orbit-normal mode
                call cross( -Sat_Coor_C, Sat_Vel_C, Rota_S2C(1:3,2))   ! Y正交于卫星位置-速度平面
                Rota_S2C(1:3,2)= Rota_S2C(1:3,2)/dsqrt(DOT_PRODUCT(Rota_S2C(1:3,2), Rota_S2C(1:3,2)))
                call cross(Rota_S2C(1:3,2), Rota_S2C(1:3,3), Rota_S2C(1:3,1))
            end if            
            ! For QZSS satellites, change to orbit normal mode phase when beta<20 degree(only QZS 1)
            ! See: http://gpsworld.com/innovation-qzs-3-and-qzs-4-join-the-quasi-zenith-satellite-system/
            if ( (System=='J') .and. (abs(beta)<20.d0) .and. (PRN==GNum+RNum+CNum+NumE+1) .and. &
             ((index(SP3File(len_trim(SP3File)-11:len_trim(SP3File)),'gbm')/=0) .or.  &
            (index(SP3File(len_trim(SP3File)-11:len_trim(SP3File)),'wum')==0)) ) then ! forQZSS, yaw-fixed, orbit-normal mode
                call cross( -Sat_Coor_C, Sat_Vel_C, Rota_S2C(1:3,2))   ! Y正交于卫星位置-速度平面
                Rota_S2C(1:3,2)= Rota_S2C(1:3,2)/dsqrt(DOT_PRODUCT(Rota_S2C(1:3,2), Rota_S2C(1:3,2)))
                call cross(Rota_S2C(1:3,2), Rota_S2C(1:3,3), Rota_S2C(1:3,1))
            end if

end subroutine