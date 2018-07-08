! ================ Bancroft ==============
! Function: 
!              This subroutine is to calculate the appropriate
!  coordinate of the station using SPP
!
! INPUTS:
!      ObsDataRaw    ObsData type
!         k                     the kth station
! Output:
!         Coor            Appropriate coordinate of the station
!
! WRITTEN BY: Yize Zhang, zhyize@163.com
! ============= End of Header =============

subroutine Bancroft(ObsDataRaw, k, Coor)
use MOD_ObsData
use MOD_constant
use MOD_GLO_Fre
use MOD_VAR
use MOD_STA
use MOD_FileID
implicit none
    type(type_ObsData) :: ObsDataRaw, ObsData
    real(8) :: Coor(3)
    integer :: Obsweek
    real(8) :: Obssec
    integer :: i, j, PRN,N, k, KK
    character(1) :: System
    real(8) :: C1,P1, P2, Range, t1, Rec_Clk, Rela, s, Sat_Coor(3), Sat_Coor0(3), Sat_Clk,toe
    real(8) :: A(40,4), L(40), V(40), dx(4), Nbb(4,4), InvN(4,4), U(4),sigma0
    real(8) :: Lat,Lon, Hgt, Sat_Vel(3)
    real(8) :: press ,tempre, rhumity, undu, MJD
    real(8) :: Lat2(2), Lon2(2), Hgt2(2), ah2(2), aw2(2),rh2(2), la2(2),undu2(2),press2(2), temp2(2),dT2(2),Tm2(2)
    logical :: If_Trop
    real(8) :: Rotation(3,3), Sat_XYZ(3), Ele, ZHD=0.d0, ZWD=0.d0, MfH, MfW, STD=0.d0
    
    ObsData=ObsDataRaw
    j=0
    If_Trop=.false.
    100 Rec_Clk=0.d0
!    Coor=0.d0
    dx=0.d0  
    200 N=0
    Nbb=0.d0;    U=0.d0;    InvN=0.d0
    A=0.d0;    L=0.d0;    dx=0.d0
    
    if (If_Trop) then
        call XYZ2BLH(Coor(1), Coor(2), Coor(3),Lat,Lon,Hgt)  ! Lat and Lon is In degree
        Rotation=reshape((/ -dsind(Lat)*dcosd(Lon),  -dsind(Lon),  dcosd(Lat)*dcosd(Lon),  &
                -dsind(Lat)*dsind(Lon), dcosd(Lon) ,dcosd(Lat)*dsind(Lon),  dcosd(Lat) ,0D0, dsind(Lat)/), (/3,3/))
        call ZTD_EGNOS(Lat, Hgt, int_doy, ZHD, ZWD)
    end if
    Obsweek=ObsData%week
    Obssec=ObsData%sow
    do i=1,ObsData%PRNS
        PRN=ObsData%PRN(i)
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
            PRN=PRN+GNum+RNum+CNum
            if (.not. SystemUsed(4)) cycle
        else if (System=="J") then  ! QZSS
            PRN=PRN+GNum+RNum+CNum+NumE
            if (.not. SystemUsed(5)) cycle
        else if (System=="I") then  ! IRNSS
            PRN=PRN+GNum+RNum+CNum+NumE+JNum
            if (.not. SystemUsed(6)) cycle
        else   ! other system
            cycle
        end if
        if (SatSelected(PRN)==0) cycle

        if ((System=="G") .or. (System=='J')) then   ! GPS or QZSS
            f1=10.23d6*154d0
            f2=10.23d6*120d0
        elseif (System=="R") then   ! GLONASS
            Kk=Fre_Chann(PRN-GNum)
            f1=(1602.0d0+Kk*0.5625d0)*1.0D6   ! f1=(1602.0d0+K*0.5625d0)*1.0d6
            f2=(1246.0d0+Kk*0.4375d0)*1.0D6
        elseif (System=="C") then   ! COMPASS
            if (freq_comb=='L1L2') then
                f1=10.23d6*152.6d0
                f2=10.23d6*118.0d0
            elseif (freq_comb=='L1L3') then
                f1=10.23d6*152.6d0
                f2=10.23d6*124.0d0
            elseif (freq_comb=='L2L3') then
                f1=10.23d6*118.0d0
                f2=10.23d6*124.0d0
            end if
        elseif (System=="E") then   ! GALILEO
            if (freq_comb=='L1L2') then   ! E1 E5a
                f1=10.23d6*154.d0
                f2=10.23d6*115.0d0
            elseif (freq_comb=='L1L3') then   ! E1 E5b
                f1=10.23d6*154.d0
                f2=10.23d6*118.d0
            elseif (freq_comb=='L2L3') then   ! E5a E5b
                f1=10.23d6*115.0d0
                f2=10.23d6*118.0d0
            end if
        elseif (System=="I") then   ! IRNSS
            f1=10.23d6*115.d0
            f2=10.23d6*243.6d0
        end if

        C1=ObsData%C1(i)
        if (freq_comb=='L1L2') then
            P1=ObsData%P1(i)
            P2=ObsData%P2(i)
        elseif (freq_comb=='L1L3') then
            P1=ObsData%P1(i)
            P2=ObsData%P3(i)
        elseif (freq_comb=='L2L3') then
            P1=ObsData%P2(i)
            P2=ObsData%P3(i)
        end if
        if ((P1 /=0.0d0) .and. (P2 /=0.0d0)) then
            Range=(f1*f1*P1-f2*f2*P2)/(f1+f2)/(f1-f2) 
        elseif (P1 /=0.0d0) then
            Range=P1
        else if (P2 /=0.0d0) then
            Range=P2
        else if (C1 /=0.0d0) then
            Range=C1
        else
            cycle
        end if
         
        ! Calculate the satellite coordinate and the distance between satellite and reciver
        t1=Range/c

        call Cal_Sat_PosClk(System, Obsweek,Obssec+ObsData%Clk_Bias-Rec_Clk,  & 
                            PRN,Coor, t1, .true., Sat_Coor0,  Sat_Coor, Sat_Vel, Sat_Clk, s, Rela,toe)                    
        if ( all(abs(Sat_Coor-9999.0d0)<0.1d0) ) cycle   ! If no data of this PRN
        if (dabs(Sat_Clk-9999.0d0)<0.1d0) cycle  ! If no this Satellite clock data

        if (If_Trop) then
            Sat_XYZ=MATMUL(Rotation,(Sat_Coor-Coor))
            Ele=dasind(Sat_XYZ(3)/dsqrt(DOT_PRODUCT(Sat_XYZ,Sat_XYZ)))
            call Map_NMF(Ele, Lat, Hgt, int_doy, MFh, MFw)
            STD=ZHD*Mfh+ZWD*Mfw
        end if

        N=N+1
        A(N,:)=0.d0
        A(N,1:4)=(/(Coor-Sat_Coor)/s,  1.d0/)
        L(N)=(Range-s-STD+Sat_Clk*c-Rec_Clk*c +Rela) 

        call Change_NEQ(Nbb, U, 4, A(N,:), L(N), "add")
    end do  ! do i=1,ObsData%PRNS

    if (N>=4) then
        call InvSqrt(Nbb, 4 ,invN)
        dx=MATMUL(invN, U)
        Coor=Coor+dx(1:3)
        Rec_Clk=Rec_Clk+dx(4)/c
        V(1:N)=MATMUL(A(1:N,:), dx)-L(1:N)
        sigma0=dsqrt(DOT_PRODUCT(V(1:N), V(1:N))/(N-4) )
        if (any(abs(dx)>1.d0)) then  ! Loop untill dx<1.d0
            goto 200
        elseif (sigma0>100.d0) then   ! if one satellite unhealth
            j=j+1
            if (j>ObsData%PRNS) then
                Coor=0.d0
                return
            end if
            ObsData=ObsDataRaw
            ObsData%C1(j)=0.d0
            ObsData%P1(j)=0.d0
            ObsData%P2(j)=0.d0
            goto 100
        elseif (.not.(If_Trop)) then
            If_Trop=.true.
            goto 200
        end if
    end if

    write(LogID,'(A15,I5,F10.1,3F15.3)') '%%%AppCoor%%%', Obsweek, Obssec,Coor

    if (all(STA%STA(k)%TrueCoor==0.d0)) then  ! Only the first epoch to calculate ZTD and compare NEU
        call XYZ2BLH(Coor(1), Coor(2), Coor(3),Lat,Lon,Hgt)  ! Lat and Lon is In degree
        STA%STA(k)%BLH(1:3)=(/Lat,Lon, Hgt/)
        STA%STA(k)%Rotation=reshape((/ -dsind(Lat)*dcosd(Lon),  -dsind(Lon),  dcosd(Lat)*dcosd(Lon),  &
                -dsind(Lat)*dsind(Lon), dcosd(Lon) ,dcosd(Lat)*dsind(Lon),  dcosd(Lat) ,0D0, dsind(Lat)/), (/3,3/))
        STA%STA(k)%Coor(1:3)=Coor+MATMUL(STA%STA(k)%NEU, STA%STA(k)%Rotation)
        STA%STA(k)%TrueCoor(1:3)=Coor
        call XYZ2BLH(STA%STA(k)%Coor(1), STA%STA(k)%Coor(2), STA%STA(k)%Coor(3),Lat,Lon,Hgt)  ! Lat and Lon is In degree
        STA%STA(k)%BLH(1:3)=(/Lat,Lon, Hgt/)
        call UTC2MJD(int_year,1,int_doy, 12, 0, 0.d0, MJD)   ! In the middle of the day
        call gpt (MJD, Lat*pi/180.d0, Lon*pi/180.d0, Hgt, press ,tempre, rhumity, undu)
        if ( (cdattype(1:3)=='NOM') .or.  (trim(cdattype)=='GPT') ) then
            call gpt (MJD, Lat*pi/180.d0, Lon*pi/180.d0, Hgt, STA%STA(k)%TROP%press ,STA%STA(k)%TROP%temp, STA%STA(k)%TROP%e, STA%STA(k)%TROP%undu)
        elseif (cdattype(1:5)=='GPT2 ') then 
            Lat2(1)=Lat*pi/180.d0; Lon2(1)=Lon*pi/180.d0; Hgt2(1)=Hgt
            call gpt2(MJD, Lat2, Lon2, Hgt2, 1, 0,press2,temp2,undu2,rh2,ah2,aw2,undu2)
        elseif (cdattype(1:7)=='GPT2_1w') then 
            Lat2(1)=Lat*pi/180.d0; Lon2(1)=Lon*pi/180.d0; Hgt2(1)=Hgt
            call gpt2_1w(MJD, Lat2, Lon2, Hgt2, 1, 0,press2,temp2,dT2,Tm2,rh2,ah2,aw2,la2,undu2)
        elseif (cdattype(1:7)=='GPT2_5w') then 
            Lat2(1)=Lat*pi/180.d0; Lon2(1)=Lon*pi/180.d0; Hgt2(1)=Hgt
            call gpt2_5w(MJD, Lat2, Lon2, Hgt2, 1, 0,press2,temp2,dT2,Tm2,rh2,ah2,aw2,la2,undu2)
        else
            write(*,*) '***ERROR : unknown ZTD input meteo data : ',cdattype
            pause
            stop
        end if
        if (cdattype(1:4)=='GPT2') then
            STA%STA(k)%TROP%press=press2(1)
            STA%STA(k)%TROP%temp=temp2(1)
            STA%STA(k)%TROP%dT=dT2(1)
            STA%STA(k)%TROP%Tm=Tm2(1)
            STA%STA(k)%TROP%e=rh2(1)
            STA%STA(k)%TROP%ah=ah2(1)
            STA%STA(k)%TROP%aw=aw2(1)
            STA%STA(k)%TROP%la=la2(1)
        end if
        ! Calculate the tropsphere ZHD and ZWD
        if (cztd(1:5)=='EGNOS') then    ! EGNOS to get ZHD and ZTD
            call ZTD_EGNOS(Lat, Hgt, int_doy, STA%STA(k)%Trop%ZHD, STA%STA(k)%Trop%ZWD)
        else if (cztd(1:4)=='SAAS') then  ! SAAS
            call ZTD_SAAS(STA%STA(k)%TROP%press ,STA%STA(k)%TROP%temp, STA%STA(k)%TROP%e, Lat*pi/180.d0, Hgt, STA%STA(k)%TROP%ZHD,STA%STA(k)%TROP%ZWD)
            if (cdattype(7:7)=='w') then
                call ZWD_Askne(STA%STA(k)%TROP%e,STA%STA(k)%TROP%Tm,STA%STA(k)%TROP%la,STA%STA(k)%TROP%ZWD)
            end if
        else if (cztd(1:5)=='UNB3 ') then  ! UNB3
            call ZTD_UNB3(Lat, Hgt, int_doy, STA%STA(k)%Trop%ZHD,STA%STA(k)%Trop%ZWD)
        elseif (cztd(1:5)=='UNB3m') then  ! UNB3m
            call ZTD_UNB3m(Lat, Hgt, int_doy, STA%STA(k)%Trop%ZHD,STA%STA(k)%Trop%ZWD)
        else
            write(*,*) '***ERROR : unknown ZTD model: ' ,cztd
            pause
            stop
        end if
    end if
    
    return
end subroutine
