! ================Cal_Rec_Clk=================
! Function:Calculate the reciver clock corrrection using the
!       presudorange. It is used to calculate the precise satellite coordinate
!
! INPUTS:
!       Obsweek        Observation week
!       Obssec           Observation seconds in a week
!       ObsData         Observation data in a specific station
!       Coor                Coordinate of the station, unit in meter
!       Rotation          Rotation matirx(3*3) from XYZ rectangular coordinate system
!                               to station fixed coordinate system
! OUTPUT:
!        Rec_Clk         receiver clock correction,unit in second
!  
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================= End of header ===============

subroutine Cal_Rec_Clk(Obsweek,Obssec,ObsData,Coor, Rotation,Rec_Clk)
use MOD_constant
use MOD_ObsData
use MOD_Glo_Fre
use MOD_FileID
use MOD_Var
use MOD_NavData
implicit none
    ! Intent in
    type(type_ObsData) :: ObsData
    integer Obsweek
    real(8) :: Obssec
    real(8) :: Coor(3)
    real(8) :: Rotation(3,3)
    ! Intent out
    real(8) :: Rec_Clk, Rec_Clk_PRN
    
    integer :: i, PRN
    real(8)  :: k, KK, Freq
    character(1) :: System
    real(8) :: P1, P2, C1, C2, Range, t1
    real(8) :: Sat_Coor(3), Sat_Coor0(3), Sat_Vel(3), s
    real(8) :: Rela, Ele
    real(8) :: Sat_Clk,toe
    
    Rec_Clk=0.0D0
    KK=0.0d0
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
            f1=f_L1
            f2=f_L2
        elseif (System=="R") then   ! GLONASS
            Freq=Fre_Chann(PRN-GNum)
            f1=(1602.0d0+Freq*0.5625d0)*1.0D6   ! f1=(1602.0d0+K*0.5625d0)*1.0d6
            f2=(1246.0d0+Freq*0.4375d0)*1.0D6
        elseif (System=="C") then   ! COMPASS
            if (freq_comb=='L1L2') then
                f1=f_B1
                f2=f_B2
            elseif (freq_comb=='L1L3') then
                f1=f_B1
                f2=f_B3
            elseif (freq_comb=='L2L3') then
                f1=f_B2
                f2=f_B3
            end if
        elseif (System=="E") then   ! GALILEO
            if (freq_comb=='L1L2') then   ! E1 E5a
                f1=f_E1
                f2=f_E5a
            elseif (freq_comb=='L1L3') then   ! E1 E5b
                f1=f_E1
                f2=f_E5b
            elseif (freq_comb=='L2L3') then   ! E5a E5b
                f1=f_E5a
                f2=f_E5b
            end if
        elseif (System=="I") then   ! IRNSS
            f1=f_L1
            f2=f_S
        end if

        C1=ObsData%C1(i)
        C2=ObsData%C2(i)
        if (freq_comb=='L1L2') then
            P1=ObsData%P1(i)
            P2=ObsData%P2(i)
            if ((P1==0.d0) .and. (C1/=0.d0)) P1=C1+DCB(PRN)*c
            if ((P2==0.d0) .and. (C2/=0.d0)) P2=C2 
            if ((P1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=P1- NavData(PRN)%Nav(2)%TGD(1)*c   ! P1(Y)
            if ((P2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=P2- (f1/f2)**2*NavData(PRN)%Nav(2)%TGD(1)*c  ! P2(Y)
            if ((P1==0.d0) .and. (C1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=C1- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL1CA*c ! L1C/A
            if ((P2==0.d0) .and. (C2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=C2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL2C*c   ! L2C
        elseif (freq_comb=='L1L3') then
            P1=ObsData%P1(i)
            P2=ObsData%P3(i)
            if ((P1==0.d0) .and. (C1/=0.d0)) P1=C1+DCB(PRN)*c
            if ((P1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=P1- NavData(PRN)%Nav(2)%TGD(1)*c   ! P1(Y)
            if ((P1==0.d0) .and. (C1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=C1- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL1CA*c ! L1C/A
            if ((P2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=P2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL5Q5*c  ! L5Q5
        elseif (freq_comb=='L2L3') then
            P1=ObsData%P2(i)
            P2=ObsData%P3(i)
            if ((P1==0.d0) .and. (C2/=0.d0)) P1=C2
            if ((P1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=P1- (f1/f2)**2*NavData(PRN)%Nav(2)%TGD(1)*c  ! P2(Y)
            if ((P1==0.d0) .and. (C2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=C2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL2C*c   ! L2C
            if ((P2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=P2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL5Q5*c  ! L5Q5
        end if
        if ((P1 /=0.0d0) .and. (P2 /=0.0d0)) then
            Range=(f1*f1*P1-f2*f2*P2)/(f1+f2)/(f1-f2)
        else
            Range=0.d0
        end if   ! if ((P1 /=0.0) .and. (P2 /=0)) then
        
        if (index(ObsCombine,"PC")==0) then
            if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
                Range=P1
            elseif ( ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) .and. (freq_comb=='L2L3')) then
                Range=P1
            else
                Range=P2
            end if
        end if
        if (Range==0.d0) cycle
            
            ! Calculate the satellite coordinate and the distance between satellite and reciver
            t1=Range/c
            call Cal_Sat_PosClk(System, Obsweek,Obssec+ObsData%Clk_Bias, PRN, Coor, t1, .true., Sat_Coor0,Sat_Coor, Sat_Vel,Sat_Clk, s, Rela,toe)
            if ( all(abs(Sat_Coor-9999.0d0)<0.1d0) ) cycle   ! If no data of this PRN
            if (dabs(Sat_Clk-9999.0d0)<0.1d0) cycle  ! If no this Satellite clock data

        Sat_Coor=MATMUL(Rotation,(Sat_Coor-Coor))
        Ele=dasind(Sat_Coor(3)/dsqrt(DOT_PRODUCT(Sat_Coor,Sat_Coor)))
        if (Ele<LimEle) cycle ! The cutoff satellite elevation
            
        KK=KK+k ! The coefficient
            ! The recvier clock corrction for each satellite
        Rec_Clk_prn=(Range - s - 5.0d0+Rela) /c +Sat_Clk  ! The trop is about 5.0 m, and rela=0.0
!            if ( (i>1) .and. (dabs(Rec_Clk - Rec_Clk_prn)>1.0D-7) .and. (dabs(Rec_Clk_prn)>1.0D-6)) then
        Rec_Clk=Rec_Clk + k/KK*(Rec_Clk_prn - Rec_Clk)
!            else ! If the diff is less than 1D-7 s or 
                !the clock correction is less than 1D-6 s, then we can ignore it.
!                exit
!            end if  ! if ((abs(Rec_Clk - R
    end do  ! do i=1,ObsData%PRNS
    
    return
end subroutine
    
    

