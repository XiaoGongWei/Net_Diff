! ============= Zero Difference ===============
! PURPOSE:
!         Form zero difference equation for each station
! INPUTS:
!      ObsData          obsdata
!          k                    the kth   station
! OUTPUT:
!         ZD                  Zero difference strcuture
!
! WRITTEN BY: Yize Zhang
! ==========================================

subroutine Zero_Difference(Obsdata, k, ZD)
use MOD_ZD
use MOD_NavHead
use MOD_ObsData
use MOD_STA
use MOD_VAR
use MOD_FileID
use MOD_CycleSlip
use MOD_Ant
use MOD_EOP
use MOD_Rotation
use MOD_GLO_Fre
use MOD_constant
use MOD_NavData
use MOD_NEQ_DP
implicit none
    type(type_ZD) :: ZD
    type(type_ObsData) :: ObsData
    integer :: k
    ! Local variables
    real(8) :: AppCoor(3),ZHD, ZWD, Lat, Lon, Hgt
    integer :: Obsweek
    real(8) :: Obssec
    real(8) :: Rec_Clk
    real(8) :: FHR, dx_SolTide(3), dx_Ocean(3), dx_CMC(3), dx_pole(3)
    real(8) :: xpm, ypm, coLat
    real(8) :: StaPCO(2)=0.d0,SatPCO(2)=0.d0,SatPCV(2)=0.d0, StaPCV(2)=0.d0, Ele_Sat

    integer :: i, j, N, PRN, PRN_S, Kk, PRNOUT(10), PRNOUTn=0, freq, N_DP
    integer(1) :: Sys, LLI1, LLI2, Slip
    character(1) :: System
    real(8) :: P1, P2, P3, C1, C2, DP1, DP2, DP3, L1, L2, L3, Range, t1,PC, LC, ZD_L1, ZD_L2, ZD_L3, EWL_amb
    real(8) :: Sat_Coor0(3),Sat_Coor(3), Sat_Vel(3), Sat_Clk, Sat_Clk_Vel, Rela, s
    real(8) :: Sat_XYZ(3), Ele, Azi, P, MJD, STD, Ion1, Ion2, Ion3, corr, dist
    real(8) :: dztd, ZTDsec
    character(100) :: line

    ZD%PRN=0;      ZD%A=0.d0
    ZD%Ele=0.d0;   ZD%P=0.d0
    ZD%P1=0.d0;    ZD%P2=0.d0
    ZD%P1CS=0.d0;    ZD%P2CS=0.d0
    ZD%L1=0.d0;    ZD%L2=0.d0
    ZD%WL=0.d0;  ZD%W4=0.d0
    ZD%EWL=0.d0; !ZD%WL_amb=99.d0


    Obsweek=ObsData%week
    Obssec=ObsData%sow
    AppCoor=STA%STA(k)%Coor(1:3)++MATMUL( STA%STA(k)%NEU, STA%STA(k)%Rotation )

    ZHD=STA%STA(k)%Trop%ZHD
    ZWD=STA%STA(k)%Trop%ZWD
    Rotation=STA%STA(k)%Rotation
    Lat=STA%STA(k)%BLH(1)
    Lon=STA%STA(k)%BLH(2)
    Hgt=STA%STA(k)%BLH(3)
    MJD=Obsweek*7.d0+44244.0d0+Obssec/86400.d0
    dist=dsqrt(DOT_PRODUCT((AppCoor-STA%STA(2)%TrueCoor),(AppCoor-STA%STA(2)%TrueCoor)))
    if ((k==2) .and. (dist>1.d3)) then ! If rover station move more than 1km, re-calculate the rotation matrix
        call XYZ2BLH(AppCoor(1), AppCoor(2), AppCoor(3),Lat,Lon,Hgt)
        Rotation=reshape((/ -dsind(Lat)*dcosd(Lon),  -dsind(Lon),  dcosd(Lat)*dcosd(Lon),  &
                    -dsind(Lat)*dsind(Lon), dcosd(Lon) ,dcosd(Lat)*dsind(Lon),  dcosd(Lat) ,0D0, dsind(Lat)/), (/3,3/))
    end if

    ! 说明(对于long baseline RTK, 这些都需要考虑)：
    ! （1）固体潮，海潮，测站PCO,PCV不需考虑，在星间差中（已改正）
    ! （2）卫星PCO,PCV（已改正）,卫星相位缠绕不需考虑，在站间差中,但如果测站太远可能有点影响？
    ! （4）测站钟差不需考虑，可以在星间差中消除，但需要估计概略值。
    ! （5）电离层公共部分可以在站间差中消除，其他受测站间距离影响
    ! （6）对流层用模型改(不改可能也可以，在站间差中消除)，其他受测站间距离影响
    
    if (If_Est_Iono) then  ! For long baseline RTK, should consider solid tide and ocean load correction
        ! 1. Calculate the solid tide correction
    !  Reference :  http://www.navipedia.net/index.php/Solid_Tides
        FHR=(ObsData%hour+ObsData%min/60.d0+ObsData%sec/3600.d0)
        dx_SolTide=0.d0
        call Tide(AppCoor,ObsData%year,ObsData%mon,ObsData%day,FHR, &
                &     MATMUL(Rota_C2T,SunCoor(1:3)), MATMUL(Rota_C2T,MoonCoor(1:3)), dx_SolTide)  ! In TRS
        
            ! 2. Calculate the ocean tide correction
        !   Reference : http://www.navipedia.net/index.php/Ocean_loading
        dx_Ocean=0.d0
        call OceanLoad('***',ObsData%year, int_doy+FHR/24.d0,STA%STA(k)%OLC%OLC, dx_Ocean) ! In NEU
        ! Orbits is always referred to geocenter, so cmc must be added to station.
        ! If l_olc_with_cmc (coefficients with CMC), then no additional CMC correction is needed.  ! Reference: A_RTK: tide_displace.f90
        ! For satellite techniques, the crust-fixed stations should include the 'geocenter motion'.
        ! For VLBI, neglect of geocenter motion have no observable consequences. ! Reference: IERS 2010,p:109
        dx_CMC=0.d0

        ! 3. Calculate the poe tide correction
        ! Displacement in east, south and radial in mm
        ! Reference : IERS Conventions (2003),pp84
        ! Reference : A_RTK MOD_others..f90
        ! Reference : http://www.navipedia.net/index.php/Pole_Tide
        dx_pole=0.d0
        xpm=MP%xp0 + MP%xp_rate*(mjd-MP%tref)/365.25d0  ! in arcseconds
        ypm=MP%yp0 + MP%yp_rate*(mjd-MP%tref)/365.25d0
        xpm=EOP%Xp*2.062648062470964d5-xpm   ! in second of arc
        ypm=ypm-EOP%Yp*2.062648062470964d5
        colat=90.d0-Lat
        dx_pole(2)=  9.d0*dcosd(colat)     *(xpm*dsind(Lon)-ypm*dcosd(Lon))     ! lambda, East
        dx_pole(1)=  - 9.d0*dcosd(2.d0*colat)*(xpm*dcosd(Lon)+ypm*dsind(Lon))   ! phi, North
        dx_pole(3)= -32.d0*dsind(2.d0*colat)*(xpm*dcosd(Lon)+ypm*dsind(Lon))   ! Up
        ! unit to meters
        dx_pole=dx_pole*1d-3

        AppCoor=AppCoor +  dx_SolTide + +MATMUL(dx_Ocean+dx_pole, Rotation) ! dx_Ocean from NEU to XYZ
    end if

    ! ********* Calculate the initial value of reciver clock corrrection *********
!    call Cal_Rec_Clk(Obsweek, Obssec, ObsData, AppCoor, Rotation, Rec_Clk)
    call Cal_Rec_Clk2(k, Obsweek, Obssec, ObsData,AppCoor, Rotation, ZHD,ZWD,Lat, Hgt, int_doy, PRNOUT, PRNOUTn, Rec_Clk)
    if (ObsData%PRNS - PRNOUTn <=3) then
        write(unit=LogID,fmt='(A5,I3,A30)') '%STA', k, 'too few satellites, skip.'
        return
    end if

    write(unit=LogID,fmt='(A5,I3,E15.6)') '%STA', k, rec_clk
    dztd=0.d0
!    read(unit=AmbID(k),fmt='(F10.4)') dztd
!    read(unit=AmbID(k),fmt='(A)') line
!    read(line(14:18),'(F5.0)') ZTDSec
!    read(line(19:25),'(F7.1)') dztd
!    if (mod(ObsSec,86400.d0)-ZTDSec<300.d0) backspace(AmbID(k))
!    dztd=dztd/1000.d0-ZHD-ZWD

    N=0
    N_DP=0
    do i=1,ObsData%PRNS
        PRN=ObsData%PRN(i)
        PRN_S=PRN
        System=ObsData%System(i)
        if (System=="G") then  !! GPS
            if (.not. SystemUsed(1)) cycle
        else if (System=="R") then  ! GLONASS
            PRN=PRN+GNum
            if (.not. SystemUsed(2)) cycle
        else if (System=="C") then  ! COMPASS
            PRN=PRN+GNum+RNum
            if (.not. SystemUsed(3)) cycle
        else if (System=="E") then  ! GALILEO
            if (.not. SystemUsed(4)) cycle
            PRN=PRN+GNum+RNum+CNum
        else if (System=="J") then  ! QZSS
            if (.not. SystemUsed(5)) cycle
            PRN=PRN+GNum+RNum+CNum+NumE
        else if (System=="I") then  ! IRNSS
            if (.not. SystemUsed(6)) cycle
            PRN=PRN+GNum+RNum+CNum+NumE+JNum
        else   ! other system
        100    cycle
        end if
        if (SatSelected(PRN)==0) cycle
        do j=1, PRNOUTn
            if (PRN==PRNOUT(j)) goto 100
        end do

        if ((System=="G") .or. (System=='J')) then   ! GPS/QZSS
            f1=f_L1
            f2=f_L2
            f3=f_L5
            Sys=1
        elseif (System=="R") then   ! GLONASS
            Kk=Fre_Chann(PRN-GNum)
            f1=(1602.0d0+Kk*0.5625d0)*1.0D6   ! f1=(1602.0d0+K*0.5625d0)*1.0d6
            f2=(1246.0d0+Kk*0.4375d0)*1.0D6
            Sys=2
        elseif (System=="C") then   ! COMPASS
            if (freq_comb=='L1L2') then
                f1=f_B1
                f2=f_B2
                f3=f_B3
            elseif (freq_comb=='L1L3') then
                f1=f_B1
                f2=f_B3
            elseif (freq_comb=='L2L3') then
                f1=f_B2
                f2=f_B3
            end if
            Sys=3
        elseif (System=="E") then   ! GALILEO
            if (freq_comb=='L1L2') then   ! E1 E5a
                f1=f_E1
                f2=f_E5
                f3=f_E5b
            elseif (freq_comb=='L1L3') then   ! E1 E5b
                f1=f_E1
                f2=f_E5b
            elseif (freq_comb=='L2L3') then   ! E5a E5b
                f1=f_E5a
                f2=f_E5b
            end if
            Sys=4
        elseif (System=="I") then   ! IRNSS
            f1=f_L1
            f2=f_S
            Sys=5
        else
            cycle
        end if
                 ! Range
        C1=ObsData%C1(i)
        C2=ObsData%C2(i)
        if (freq_comb=='L1L2') then
            P1=ObsData%P1(i)
            P2=ObsData%P2(i)
            P3=ObsData%P3(i)
            if ((P1==0.d0) .and. (C1/=0.d0)) P1=C1+DCB(PRN)*c
!            if ( (C2/=0.d0)) P2=C2 
            if ((P2==0.d0) .and. (C2/=0.d0)) P2=C2 
            if ((P1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=P1- NavData(PRN)%Nav(2)%TGD(1)*c   ! P1(Y)
            if ((P2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=P2- (f1/f2)**2*NavData(PRN)%Nav(2)%TGD(1)*c  ! P2(Y)
            if ((P1==0.d0) .and. (C1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=C1- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL1CA*c ! L1C/A
            if ((P2==0.d0) .and. (C2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=C2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL2C*c   ! L2C
            L1=ObsData%L1(i)
            L2=ObsData%L2(i)
!            if (k==2 .and. PRN>GNum+RNum .and. PRN<GNum+RNum+6) then
!                L1=L1-0.5d0
!                L2=L2-0.5d0
!            end if
            L3=ObsData%L3(i)
            LLI1=ObsData%LLI1(i)
            LLI2=ObsData%LLI2(i)
            if (L2==0.d0 .and. ObsData%L2C(i)/=0.d0) then
                L2=ObsData%L2C(i)
                LLI2=ObsData%LLI2C(i)
            end if
            DP1=ObsData%D1(i)
            DP2=ObsData%D2(i)
        elseif (freq_comb=='L1L3') then
            P1=ObsData%P1(i)
            P2=ObsData%P3(i)
            if ((P1==0.d0) .and. (C1/=0.d0)) P1=C1+DCB(PRN)*c
            if ((P1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=P1- NavData(PRN)%Nav(2)%TGD(1)*c   ! P1(Y)
            if ((P1==0.d0) .and. (C1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=C1- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL1CA*c ! L1C/A
            if ((P2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=P2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL5Q5*c  ! L5Q5
            L1=ObsData%L1(i)
            L2=ObsData%L3(i)
            LLI1=ObsData%LLI1(i)
            LLI2=ObsData%LLI3(i)
            DP1=ObsData%D1(i)
            DP2=ObsData%D3(i)
        elseif (freq_comb=='L2L3') then
            P1=ObsData%P2(i)
            P2=ObsData%P3(i)
            if ((P1==0.d0) .and. (C2/=0.d0)) P1=C2
            if ((P1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=P1- (f1/f2)**2*NavData(PRN)%Nav(2)%TGD(1)*c  ! P2(Y)
            if ((P1==0.d0) .and. (C2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=C2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL2C*c   ! L2C
            if ((P2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=P2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL5Q5*c  ! L5Q5
            L1=ObsData%L2(i)
            L2=ObsData%L3(i)
            LLI1=ObsData%LLI2(i)
            LLI2=ObsData%LLI3(i)
            DP1=ObsData%D2(i)
            DP2=ObsData%D3(i)
            if (L1==0.d0 .and. ObsData%L2C(i)/=0.d0) then
                L1=ObsData%L2C(i)
                LLI1=ObsData%LLI2C(i)
            end if
        end if
        if ((P1 /=0.0d0) .and. (P2 /=0.0d0)) then
            Range=(f1*f1*P1-f2*f2*P2)/(f1+f2)/(f1-f2) ! Ionospheric-free combination
        else if (P1 /=0.0d0) then
            Range=P1
        else if (C1 /=0.0d0) then
            P1=C1
            Range=P1
        else if (P2 /=0.0d0) then
            Range=P2
        else
            cycle
        end if   ! if ((P1 /=0.0) .and. (P2 /=0))
        t1=Range/c
        
        ! ********Satellite coordinate and clcok ********
        call Cal_Sat_PosClk(System, Obsweek, Obssec+ObsData%Clk_Bias-Rec_Clk,PRN, AppCoor, t1, .true., Sat_Coor, Sat_Vel, Sat_Clk, s, Rela)
        if ( all(dabs(Sat_Coor-9999.0d0)<0.1d0) ) cycle   ! If no data of this PRN
        if (dabs(Sat_Clk-9999.d0)<0.1d0 ) cycle
               
        Sat_XYZ=MATMUL(Rotation,(Sat_Coor-AppCoor))
        Ele=dasind(Sat_XYZ(3)/dsqrt(DOT_PRODUCT(Sat_XYZ,Sat_XYZ)))
        Azi=datand(Sat_XYZ(2)/Sat_XYZ(1))
        if (Sat_XYZ(1)<0.d0) then
            Azi=Azi+180.d0
        elseif (Sat_XYZ(2)<0.d0) then
            Azi=Azi+360.d0
        end if

        if (Ele<LimEle) then
            write(LogID,'(A6,1X,A1,I2,F8.2,A20)')  'PRN', System, PRN_S, Ele, 'elevation too low'
            cycle
        elseif (Ele>30.d0) then
            P=1.d0/(0.5d0+0.5d0/sind(Ele))**2 !1.d0   ! 0.3 is refer from RTKLIB
        elseif (Ele>LimEle) then
            P=1.d0/(0.5d0+0.5d0/sind(Ele))**2 !4*dsind(Ele)**2   !
        end if

        ! **************** Satellite and Receiver PCO PCV ***************
        StaPCO(1)=DOT_PRODUCT(Ant(GNum+RNum+CNum+NumE+JNum+INum+k)%PCO(1:3,1), (/cosd(Ele)*cosd(Azi), cosd(Ele)*sind(Azi), sind(Ele)/) ) ! in L1
        StaPCO(2)=DOT_PRODUCT(Ant(GNum+RNum+CNum+NumE+JNum+INum+k)%PCO(1:3,2), (/cosd(Ele)*cosd(Azi), cosd(Ele)*sind(Azi), sind(Ele)/) )  ! in L2
        StaPCV=0.d0; SatPCV=0.d0; SatPCO=0.d0;
        if (allocated(Ant(GNum+RNum+CNum+NumE+JNum+INum+k)%PCV)) then
            call PCV_Corr(GNum+RNum+CNum+NumE+JNum+INum+k, 90.d0-Ele, Azi, StaPCV)  ! Station PCV, zenith and azimuth depended
        end if
         call S2C(SunCoor(1:3), Sat_Coor, Sat_Vel, System, PRN)
        if (Orbit=="SP3") then
             Ele_Sat=dasind(sind(Ele+90.d0)*6378137.d0/dsqrt(DOT_PRODUCT(Sat_Coor,Sat_Coor)))
!            SatPCO(1)=Ant(PRN)%PCO(3,1)*cosd(Ele_Sat)  ! simply consider the Z-PCO difference, incorrect!!!!
!            SatPCO(2)=Ant(PRN)%PCO(3,2)*cosd(Ele_Sat)
            Sat_Coor0=Sat_Coor+MATMUL(Rota_C2T, MATMUL(Rota_S2C,Ant(PRN)%PCO(1:3,1)))
            SatPCO(1)=dsqrt(DOT_PRODUCT((Sat_Coor-AppCoor),(Sat_Coor-AppCoor)))-dsqrt(DOT_PRODUCT((Sat_Coor0-AppCoor),(Sat_Coor0-AppCoor)))
            Sat_Coor0=Sat_Coor+MATMUL(Rota_C2T, MATMUL(Rota_S2C,Ant(PRN)%PCO(1:3,2)))
            SatPCO(2)=dsqrt(DOT_PRODUCT((Sat_Coor-AppCoor),(Sat_Coor-AppCoor)))-dsqrt(DOT_PRODUCT((Sat_Coor0-AppCoor),(Sat_Coor0-AppCoor)))
            call PCV_Corr(PRN, Ele_Sat, 0.d0, SatPCV)  ! Satllite PCV, zenith depended
        end if
        ! Note:
        ! The difference of satellite satellite zenith angle on double difference is dis(km)/26000*180/pi
        ! -------For 1000km baseline, the satellite zenith angle difference is 2.2deg
        ! For PCO of 2m(precise orbit), it is 7cm; For PCV, it is within 1cm
        ! -------For 100km baseline, the satellite zenith angle difference is 0.22deg
        ! For PCO of 2m(precise orbit), it is 0.7cm; For PCV, it can be ignored
        ! ----------------------------------------------------------------------------------------------------

        ! ************STD(Slant Tropsphere Delay)*************
        if  (cmap(1:4)=='SAAS') then
            STA%STA(k)%Trop%map_dry=1.d0/dsind(Ele)
            STA%STA(k)%Trop%map_dry=STA%STA(k)%Trop%map_wet
        elseif (cmap(1:3)=='NMF') then
            call Map_NMF(Ele, Lat, Hgt, int_doy, STA%STA(k)%Trop%map_dry, STA%STA(k)%Trop%map_wet)
        elseif (cmap(1:3)=='GMF') then
            !pause
            call Map_GMF(MJD,Lat*pi/180.d0,Lon*pi/180.d0, Hgt,  (90.d0-Ele)*pi/180.d0, STA%STA(k)%Trop%map_dry, STA%STA(k)%Trop%map_wet)
        elseif ( (cmap(1:4)=='VMF1') .and. (cdattype(1:4)=='GPT2') )  then
            call vmf1_ht(STA%STA(k)%Trop%ah,STA%STA(k)%Trop%aw, MJD, Lat*pi/180.d0, Hgt, (90.d0-Ele)*pi/180.d0, STA%STA(k)%Trop%map_dry, STA%STA(k)%Trop%map_wet)
            ! call vmf1(STA%STA(k)%Trop%ah,STA%STA(k)%Trop%aw, MJD, Lat*pi/180.d0, (90.d0-Ele)*pi/180.d0, STA%STA(k)%Trop%map_dry, STA%STA(k)%Trop%map_wet)
        else
            write(*,*) '***ERROR : unknown ZTD mapping  function or wrong combination: ',cmap,cdattype
            pause
            stop
        end if
        STD=STA%STA(k)%Trop%map_wet*(ZWD+dztd)+STA%STA(k)%Trop%map_dry*ZHD
!        call UNB3(Lat, Ele,  Hgt, int_doy, STD)

        if (iontype==1) then
            if (all(NavHead%Alpha/=0.d0)) then
                call Klobuchar(Lat, Lon, Ele, Azi, Obssec, NavHead%Alpha,NavHead%Beta,Ion1)
            else
                Ion1=0.d0
            end if
        elseif (iontype==2) then
            call GIM(Lat, Lon, Ele, Azi, ObsWeek, ObsSec,Ion1)
        end if
        Ion2=Ion1*f1**2/f2**2
        Ion3=Ion1*f1**2/f3**2

        ! ============= Doppler Velocity ==============
        if (Combination(3) .and. k==2 .and. Dp1/=0.d0) then
            N_DP=N_DP+1
            Sat_Vel=Sat_Vel + dsin(omg)*(/Sat_Coor(2), -Sat_Coor(1), 0.d0/) ! account for earth rotation velocity,  in ecef
            NEQ_DP%PRN(N_DP)=PRN
            NEQ_DP%P(N_DP)=P/sigDP
            NEQ_DP%A(N_DP, :) = (/-(AppCoor-Sat_Coor)/s, -1.d0/)*P/sigDP
            Sat_Clk_Vel=0.d0 ! Sat clock velocity is usually less than 1e-10s/s, so the error is less than 0.03m/s, this can be ignored
            NEQ_DP%L(N_DP)=(Dp1*c/f1+Sat_Clk_Vel*c-DOT_PRODUCT((AppCoor-Sat_Coor)/s,Sat_Vel) )*P/sigDP
        end if

        if ( (L1==0.d0) .and. (L2==0.d0) ) then
            cycle
        elseif ( (L1==0.d0) .or. (L2==0.d0) ) then
            if ( (a1*f1+a2*f2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) .or. If_Est_Iono ) then   ! If dual-frequency combination
!                CycleSlip(k)%Slip(PRN)=1
                cycle   ! better than cycleslip=1
            elseif ( (a1*f1+a2*f2/=0.d0) .and. (L1==0.d0) ) then   ! If single-frequency combination
                cycle
            elseif ( (b1*f1+b2*f2/=0.d0) .and. (L2==0.d0) ) then   ! If single-frequency combination
                cycle
            end if
        end if

        ! If not instantaneous AR
        call Cycle_Slip_Detect(k, Obsweek, Obssec, P1, P2, L1, L2, LLI1, LLI2, PRN, Ele, Slip)
        if (Var_smooth=='y' .or. Var_smooth=='Y') then ! If smooth pseudorange
            ! Smooth the pseudo-range, this is very important for instanteous AR
            ! It is better to use after station difference, for it eliminates the ionosphere effect and can easily applied in single frequency
            if ((index(Smooth_Method,"Hatch") /=0)) then
                call Hatch_Filter(PRN, P1, P2, Range, L1*c/f1, L2*c/f2, Ele, k, Slip)
            elseif ((index(Smooth_Method,"Dop") /=0)) then
                if (k==1) then
                    call Hatch_Filter(PRN, P1, P2, Range, L1*c/f1, L2*c/f2, Ele, k, Slip)
                else
                    call Doppler_Filter(Obssec, PRN, P1, P2, Range, DP1, DP2,k)
                end if
            end if
        end if
        if (Slip==1) CycleSlip(k)%Slip(PRN)=1  ! Record the slip flag, this will include cycle slip information of previous information
        
!        if (k==2 .and.  sys==3) then ! BeiDou
!            P1=P1+0.68d0
!            P2=P2+0.86d0
!            L1=L1-2.44d0/c*f1
!            L2=L2-2.66d0/c*f2
!        elseif (k==2 .and.  sys==4) then ! Galileo
!            P1=P1-0.88d0
!            P2=P2+1.19d0
!            L1=L1-3.04d0/c*f1
!            L2=L2-3.26d0/c*f2
!        end if

        ! ****** Form zero difference ********
        N                    =    N+1
        ZD%Sys(N)   =    Sys
        ZD%System(N)   =    System
        ZD%PRN(N)   =    PRN
        ZD%PRN_S(N)   =    PRN_S
        ZD%Ele(N)     =    Ele
        ZD%P(N)        =    P
        ZD%A(N,1:3)  =    (/(AppCoor-Sat_Coor)/s /)
        if (TropLen/=0.d0 .and. If_Est_Iono) ZD%A(N, 4) = STA%STA(k)%Trop%map_wet ! Only when estimate ionosphere, usually for long baseline
        if (System=='R' .and. GloParaNum>0) then
            ZD%A(N, ParaNum-GloParaNum+Kk+8)=1.d0
        end if
        corr                =    s+STD-Sat_Clk*c+Rec_Clk*c-Rela 
        ZD%Corr(N)   =        STD-Sat_Clk*c+Rec_Clk*c-Rela
        ZD%s(N)        =     s
        if (P1/=0.d0)      ZD%P1(N)=P1-corr-Ion1 + StaPCO(1) - StaPCV(1)+ SatPCO(1) - SatPCV(1)
        if (P2/=0.d0)      ZD%P2(N)=P2-corr-Ion2 + StaPCO(2) - StaPCV(2)+ SatPCO(2) - SatPCV(2)
        ZD%P1CS(N)=ZD%P1(N);    ZD%P2CS(N)=ZD%P2(N)

        if ( (P1/=0.d0) .and. (P2/=0.d0)  .and. (b1*f1+b2*f2==0.d0) .and. (mod(a1,1.d0)/=0.d0) ) then   ! PC combination
            ZD%P1(N)=(f1**2*ZD%P1(N)-f2**2*ZD%P2(N))/(f1**2  -f2**2)
            ZD%P2(N)=0.d0

!            PC=(f1**2*ZD%P1(N)-f2**2*ZD%P2(N))/(f1**2  -f2**2)
!            if (P3/=0.d0)   ZD%P2(N)=P3-corr-Ion1
!            ZD%P2(N)=(f1**2*ZD%P1(N)-f3**2*ZD%P2(N))/(f1**2  -f3**2)
!            ZD%P1(N)=PC
        elseif  (a2==0.d0 .and. b2==0.d0) then   ! only P1
            ZD%P2(N)=0.d0
        elseif  (a1==0.d0 .and. b1==0.d0) then   ! only P2
            ZD%P1(N)=0.d0
        end if

!        EWL_amb=L1 - L2 - (P1*f1 + P2*f2)/(f1+f2)*(f1-f2)/c ! L1, L2 here is in cycle
        if (If_Est_Iono .and. P2/=0.d0 .and. P3/=0.d0 .and. L2/=0.d0 .and. L3/=0.d0) then
            EWL_amb=L2 - L3 - (P2*f2 + P3*f3)/(f2+f3)*(f2-f3)/c ! L2, L3 here is in cycle
            ZD%EWL_amb(N)=EWL_amb
            ZD%EWL(N)=(c*L2+Ion2*f2-c*L3-Ion3*f3)/(f2 - f3)-corr
        end if

        if (CycleSlip(k)%Slip(PRN)==1) then
            ZD%amb0(PRN,1)=real(nint(L1-Range/c*f1))
            ZD%amb0(PRN,2)=real(nint(L2-Range/c*f2))
            ZD%amb0(PRN,3)=real(nint(L3-Range/c*f3))
!            if (k==2 .and.  sys==1 .and. PRN==17) then ! GPS
!                ZD%amb0(PRN,1)=ZD%amb0(PRN,1)+4.d0
!                ZD%amb0(PRN,2)=ZD%amb0(PRN,2)+3.d0
!            end if
        elseif (IF_TC) then
        end if
            ZD%amb1(PRN,1)=real(nint(L1-Range/c*f1))  ! For tightly combined multi-system RTK
            ZD%amb1(PRN,2)=real(nint(L2-Range/c*f2))
            ZD%amb1(PRN,3)=real(nint(L3-Range/c*f3))
        if (L1/=0.d0)      ZD_L1=(L1-ZD%amb0(PRN,1))*c/f1-corr+Ion1 + StaPCO(1) - StaPCV(1)+ SatPCO(1) - SatPCV(1)    ! in distance
        if (L2/=0.d0)      ZD_L2=(L2-ZD%amb0(PRN,2))*c/f2-corr+Ion2 + StaPCO(2) - StaPCV(2)+ SatPCO(2) - SatPCV(2)    ! in distance
        if (L3/=0.d0)      ZD_L3=(L3-ZD%amb0(PRN,3))*c/f3-corr+Ion3 + StaPCO(2) - StaPCV(2)+ SatPCO(2) - SatPCV(2)    ! in distance
        if ( (a1/=0.d0) .and. (a2/=0.d0) ) then  ! If two frequency combination
            if ( (L1/=0.d0) .and. (L2/=0.d0) ) then
                ZD%WL(N)=(a1*f1*ZD_L1+a2*f2*ZD_L2)/(a1*f1+a2*f2)
!                if (CycleSlip(k)%CS(PRN)%arcLengthMW>=5) then  ! When accumulate 5 epoches, record Wide Lane ambiguity just for rounding
                if ( CycleSlip(k)%Slip(PRN)/=1 .and. Ele>LimEle) then
                      ! Just for test, not very good, because of the wrong rounding integer from code multipath
!                    ZD%WL_amb(PRN)=CycleSlip(k)%CS(PRN)%nMWmean
                    ZD%WL_amb_n(PRN)= ZD%WL_amb_n(PRN)+1.d0 ! Wide Lane ambiguity, in cycle, includes DCB
                    ZD%WL_amb(PRN)=(L1 - L2 - (P1*f1 + P2*f2)/(f1+f2)*(f1-f2)/c)/ZD%WL_amb_n(PRN) + (1.d0-1.d0/ ZD%WL_amb_n(PRN))*ZD%WL_amb(PRN)
                else  ! if new cycle slip
                    ZD%WL_amb(PRN)=99.d0
                    ZD%WL_amb_n(PRN)=0.d0
                end if
                if ((b1==0.d0) .and. (b2==0.d0)) then  ! If only WL combination
                    ZD%L1(N)=(a1*f1*ZD_L1+a2*f2*ZD_L2)/(a1*f1+a2*f2)
                elseif (b1==0.d0) then
                    ZD%L1(N)=ZD_L1
                elseif (b2==0.d0) then
                    ZD%L1(N)=ZD_L2
                end if
            end if
        elseif ( (L1/=0.d0) .and. (a1/=0.d0) ) then  ! If L1 frequency combination
            ZD%WL(N)=ZD_L1
            ZD%L1(N)=ZD_L1
        elseif ( (L2/=0.d0) .and. (a2/=0.d0) ) then  ! If L2 frequency combination
            ZD%WL(N)=ZD_L2
            ZD%L1(N)=ZD_L2
        end if
        if (If_Est_Iono .and. IonoNum>0) then
            if (L1/=0.d0) then
                ZD%L1(N)=ZD_L1
            end if
            if (L2/=0.d0) then
                ZD%L2(N)=ZD_L2
            end if
            if ( (L1/=0.d0) .and. (L3/=0.d0) ) then
!                ZD%W4(N)=(f1*ZD_L1 - f3*ZD_L3)/(f1 - f3)  ! If triple frequency ambiguity, then this can be can be applied 
            end if
        end if
        if ( (b1/=0.d0) .and. (b2/=0.d0) ) then  ! If two frequency combination
            if ( (L1/=0.d0) .and. (L2/=0.d0) ) then
                ZD%W4(N)=(b1*f1*ZD_L1+b2*f2*ZD_L2)/(b1*f1+b2*f2)
                if ((a1==0.d0) .and. (a2==0.d0)) then  ! If only W4 combination
                    ZD%L2(N)=(b1*f1*ZD_L1+b2*f2*ZD_L2)/(b1*f1+b2*f2)
                elseif (a1==0.d0) then
                    ZD%L2(N)=ZD_L1
                elseif (a2==0.d0) then
                    ZD%L2(N)=ZD_L2  
                end if
            end if
        elseif ( (L1/=0.d0) .and. (b1/=0.d0) ) then  ! If L1 frequency combination
            ZD%W4(N)=ZD_L1
            ZD%L2(N)=ZD_L1
        elseif ( (L2/=0.d0) .and. (b2/=0.d0) ) then  ! If L2 frequency combination
            ZD%W4(N)=ZD_L2
            ZD%L2(N)=ZD_L2
        end if
!        if ( (L1/=0.d0) .and. (L2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) ) ZD%W4(N)=(b1*f1*ZD%L1(N)+b2*f2*ZD%L2(N))/(b1*f1+b2*f2)
!        if ( (L1/=0.d0) .and. (L3/=0.d0) .and. (a1*f1+a2*f2/=0.d0)  ) ZD%W4(N)=(a1*f1*ZD%L1(N)+a2*f2*LC)/(a1*f1+a2*f2)
        write(unit=LogID,fmt='(A6,1X,A1,I2,2F8.2,2F13.2,E15.7,3F7.2, 2X, 2F7.2,4F13.2)') 'PRN',System,PRN_S,Ele, Azi,Range,  &
                             s  , Sat_Clk , STD, Ion1, rela, ZD%P1(N), ZD%P2(N),ZD%L1(N), ZD%L2(N),ZD%amb0(PRN,1), ZD%amb0(PRN,2) ! ZD%WL(N), ZD%W4(N)
     end do
     ZD%PRNS =  N
     ZD%week =  ObsData%week
     ZD%sow   =  ObsData%sow
     if (k==2) then
        NEQ_DP%PRNS=N_DP
     end if

     return
end subroutine