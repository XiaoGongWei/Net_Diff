! ====================  Phase_windup  ===================
!
! PURPOSE:
!        Phase windup correction.
! 
!  REFERENCE:
!      http://www.navipedia.net/index.php/Carrier_Phase_Wind-up_Effect
!
! INPUTS:
!      Rec_Coor      staiton receiver coordinate, unit in meter
!      Sat_Coor       satellite receiver coordinate, unit in meter
!      iamb             cycle slip flag
!                                 iamb=0: no cycle slip
!                                 iamb=1: cycle slip
!      dx_windup_previous    pahse windup in the previous epoch, unit in cycle
!
! OUTPUTS:
!      dx_windup      pahse windup, unit in cycle
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ===================== End of header ======================

subroutine Phase_windup(Rec_Coor, Sat_Coor, iamb,dx_windup, dx_windup_previous)
use MOD_Rotation
implicit none
    ! Intent in
    real(8) :: Rec_Coor(3), Sat_Coor(3)
    integer(1) :: iamb
    ! Intent out
    real(8) :: dx_windup, dx_windup_previous
    ! Local variables
    real(8) :: Rec_Coor_CRS(3), Sat_Coor_CRS(3)
    real(8) ::  unit_sat2rec(3)
    real(8) :: tr(3),x_rec(3), y_rec(3), x_sat(3), y_sat(3)
    real(8) :: CrossRec(3), CrossSat(3)
    real(8) :: dRec(3), dSat(3)
    real(8) :: xtmp
    real(8), parameter :: pi=3.1415926535897932d0
    real(8) :: sense(3)
    
    x_rec=Rotation(1, 1:3)   ! unit vector from NEU to TRS(ECEF), x=north
    y_rec= -Rotation(2,1:3)  ! y=west
    x_rec=Rotation(2,1:3)   ! x=east, 也可以，但是差了0.25周，吸收到钟差里
    y_rec=Rotation(1, 1:3)  ! y=north
    x_rec=matmul(Rota_T2C,x_rec)   ! unit vector of station antenna in CRS
    y_rec=matmul(Rota_T2C,y_rec)  
    x_sat=Rota_S2C(1:3,1)   ! unit vector of satellite antenna in CRS
    y_sat=Rota_S2C(1:3,2)
    
    Rec_Coor_CRS=matmul(Rota_T2C, Rec_Coor) ! Station coordinate in CRS
    Sat_Coor_CRS=MATMUL(Rota_T2C, Sat_Coor)  ! Satellite coordinate in CRS
    unit_sat2rec=(Rec_Coor_CRS - Sat_Coor_CRS)/dsqrt(DOT_PRODUCT((Rec_Coor_CRS-Sat_Coor_CRS),(Rec_Coor_CRS-Sat_Coor_CRS)) )
    tr=unit_sat2rec
    call cross(tr, y_rec, CrossRec)
    call cross(tr, y_sat, CrossSat)
    
    dRec=x_rec - tr*DOT_PRODUCT(tr,x_rec) + CrossRec
    dSat=x_sat - tr*DOT_PRODUCT(tr,x_sat) - CrossSat  ! 第二项是加还是减？？不过影响在毫米级
    
    ! phase windup in cycles
    xtmp=DOT_PRODUCT(dRec,dSat)/dsqrt(DOT_PRODUCT(dRec,dRec)*DOT_PRODUCT(dSat,dSat))
    if (xtmp>1.d0) then    ! 这是什么意思？？？？
        xtmp=1.d0
    else if(xtmp<-1.d0) then
        xtmp=-1.d0
    end if
    
    dx_windup=acos(xtmp)/2.d0/pi  ! In circles
    call cross(dSat, dRec, sense)
    if (DOT_PRODUCT(tr,sense)<0.d0)  dx_windup= - dx_windup
    
    if (iamb==1) dx_windup_previous=0.d0
    ! Add the previous value
    dx_windup=nint(dx_windup_previous - dx_windup)+dx_windup

    dx_windup_previous=dx_windup
    
    return
end subroutine