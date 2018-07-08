!  ================  XYZ2BLH  ============
!
! PURPOSE:
!    Transfer from XYZ to BLH
! 
! REFERENCES:
!    http://www.navipedia.net/index.php/Ellipsoidal_and_Cartesian_Coordinates_Conversion
! 
! INPUTS:
!       X                   coordinate of X, uniit in meter
!       Y                   coordinate of Y, uniit in meter
!       Z                   coordinate of Z, uniit in meter
!
! OUTPUTS:
!       B                  latitude,unit in degree
!       L                  longitude,unit in degree
!       H                  height,unit in meter
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!===============  End of Header  ============

subroutine XYZ2BLH(X,Y,Z,B,L,H)
use MOD_constant
implicit none
!! Elements of WGS84 ellipsolide
!real(8),parameter :: axis=6378137d0  ! major semi-axis(meter), GGSP standard
!real(8), parameter:: ecc2=0.0066943800229d0  ! square of eccentricity
!real(8), parameter :: pi=3.1415926535897932d0    ! pi
! Error tolerance
real(8), parameter :: eB=1.0D-12, eH=1.0D-5

    real(8) :: X, Y, Z
    real(8):: B,L,H
    real(8) :: p, dH, dB, B0, H0,v
    
    dH=1.0d0
    dB=1.0d0

    p=dsqrt(X*X+Y*Y)
    B=datan2(Z, p/(1.0d0 - ecc2))
    do while ((dB>eB) .or. (dH>eH))
        B0=B
        H0=H
        v=axis/ dsqrt(1.0d0- ecc2*dsin(B)*dsin(B))
        H=p*dcos(B)+Z*dsin(B)-(axis*axis)/v
        B=datan2(Z, p*(1.0d0-ecc2*v/(v+H)))
        dB=abs(B-B0)
        dH=abs(H-H0)
    end do
    L=datan2d(Y,X)
    if (L<0.d0) L=L+360.d0
    B=B*180.d0/pi
 
end subroutine