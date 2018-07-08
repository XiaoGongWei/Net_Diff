! ================== GIM ===================
! PURPOSE:
!         This function is to compute the ionospheric 
!       delay by using the IONEX file
! 
! INPUTS:
!      Lat           lattitude, in deg
!      Lon          longitude, in deg
!      Ele           elevation, in deg
!      Azi           azimuth, in deg
!      Week       week
!      Sow          second of week

! OUTPUT:
!      ION_delay             ION_delay

! Written by: Yize Zhang
!============== End of Header =============

subroutine  GIM(Lat, Lon, Ele, Azi, Week, Sow, ION_delay)
use MOD_IONHead
use MOD_IONData
!use MOD_constant
implicit none
    real(8) :: Lat ,Lon ,Ele, Azi
    integer :: week
    real(8) :: sow
    real(8) :: ION_delay
    ! Local variables
    real(8) :: rp, ap
    real(8) :: lat_i, lon_i, lon_i_init, EA
    real(8) :: lat_g, lon_g
    integer :: i, j, k, n
    real(8) :: dT, dLat, dLon
    real(8) :: IRI_T_Lat_Lon(2),IRI_T_Lat(2), IRI_T(2)
    integer :: iri_i, iri_j, iri_k
    real(8) :: sinz, cosz

    ! write(*,*) dsin(1.d0), sin(1.d0), sind(90.d0), asind(1.d0)
    lat_i=Lat-11.54544*sind(2*Lat)
    lon_i=Lon

!    EA=Ele/180.d0   ! elevation angle in semicircles
!    EA=0.0137d0/(EA+0.11d0)-0.022d0      !  Earth Centered angle
!
!    lat_i=Lat/180.d0+EA*cosd(Azi)    !  Subionospheric lat, in semicircles.
!    if (lat_i>0.416d0) then
!        lat_i=0.416d0
!    elseif (lat_i<-0.416d0) then
!        lat_i=-0.416d0
!    end if
!
!    lon_i=Lon/180.d0+EA*sind(Azi)/cos(lat_i*3.1415926d0)   ! Subionospheric long, in semicircles.
!    lat_i=lat_i*180.d0
!    lon_i=lon_i*180.d0

    ! pierce point position  ! referecne from ionppp in RTK_LIB
    rp=6371.d0/6821.d0*cosd(Ele)  !  正弦定理，穿刺点与测站和地心的天底角的正弦
    ap=90-Ele-asind(rp)           ! in degree    穿刺点与测站的地心角
    lat_i=asind( sind(Lat)*cosd(ap)+cosd(Lat)*sind(ap)*cosd(Azi) )
    lon_i_init=Lon+asind( sind(ap)*sind(Azi)/cosd(Lat_i) ) 
    
    iri_i=0
    do i=1, IONHead%Maps
        dT=(week-IONData(i)%week)*604800.d0+(sow-IONData(i)%sow)
        if (abs(dT)<=IONHead%INTERVAL) then
            iri_i=iri_i+1
            iri_j=0
            ! Rotate maps adding latitude proportional to time difference between maps 
            lon_i=lon_i_init+dT/86400.d0*360.d0
            if (lon_i>180.d0) lon_i=lon_i-360.d0      ! Longitude transfer to -180~180
            do j=1, IONHead%nLat
                dLat=IONHead%Lat1+(j-1)*IONHead%dLat-Lat_i
                if (abs(dLat)<=abs(IONHead%dLat)) then
                    iri_j=iri_j+1
                    iri_k=0
                    do k=1,  IONHead%nLon
                        dLon=IONHead%Lon1+(k-1)*IONHead%dLon-Lon_i
                        if (abs(dLon)<=abs(IONHead%dLon)) then
                            iri_k=iri_k+1
                            IRI_T_Lat_Lon(iri_k)=IONData(i)%TEC(j,k)
                        end if 
                        if (iri_k==2) then
                            IRI_T_Lat(iri_j)=IRI_T_Lat_Lon(2)+abs(dLon)/abs(IONHead%dLon)*(IRI_T_Lat_Lon(1)-IRI_T_Lat_Lon(2))
                            exit
                        end if
                    end do
                end if
                if (iri_j==2) then
                    IRI_T(iri_i)=IRI_T_Lat(2)+abs(dLat)/abs(IONHead%dLat)*(IRI_T_Lat(1)-IRI_T_Lat(2))
                    exit
                end if
            end do
        end if
        if (iri_i==2) then
            ION_delay=IRI_T(2)+abs(dT)/IONHead%INTERVAL*(IRI_T(1)-IRI_T(2))
            exit
        end if
    end do

    ! ???????
    ION_delay=ION_delay*10.d0**(IONHead%EXPONENT)*0.4028d0/1.57542d0/1.57542d0  ! See www.navipedia.net/index.php/NeQuick_Ionospheric_Model
    sinz=6371.d0/(6371.d0+IonHead%Hgt1)*sind(0.9782d0*(90.d0-ele)) ! 正弦定理，why 0.9782?
    sinz=6371.d0/(6371.d0+IonHead%Hgt1)*sind((90.d0-ele)) ! 正弦定理
    cosz=dsqrt(1-sinz*sinz)   ! From RTK_LIB ionex.c iondelay
    ION_delay=ION_delay/cosz

    return
end subroutine