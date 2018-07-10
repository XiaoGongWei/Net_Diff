!  ================  BLH2XYZ  ============
!
! PURPOSE:
!    Transfer from BLH to XYZ
! 
! REFERENCES:
!    http://www.navipedia.net/index.php/Ellipsoidal_and_Cartesian_Coordinates_Conversion
! 
! INPUTS:
!       B                  latitude,unit in degree
!       L                  longitude,unit in degree
!       H                  height,unit in meter
!
! OUTPUTS:
!       X                   coordinate of X, uniit in meter
!       Y                   coordinate of Y, uniit in meter
!       Z                   coordinate of Z, uniit in meter
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!===============  End of Header  ============

subroutine BLH2XYZ(B,L,H,X,Y,Z)
use MOD_constant
implicit none
!! Elements of WGS84 ellipsolide
!real(8),parameter :: axis=6378137d0  ! major semi-axis(meter), GGSP standard
!real(8), parameter:: ecc2=0.0066943800229d0  ! square of eccentricity
!real(8), parameter :: pi=3.1415926535897932d0    ! pi

    real(8):: B,L,H
    real(8) :: X, Y, Z
    ! Local variables
    real(8) :: N

    N=axis/sqrt(1-ecc2*(sind(B)**2))
    X=(N+H)*cosd(B)*cosd(L)
    Y=(N+H)*cosd(B)*sind(L)
    Z=(N*(1-ecc2)+H)*sind(B)

    return
end

