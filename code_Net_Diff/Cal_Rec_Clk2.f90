! ================ Cal_Rec_Clk2 =================
!
! PURPOSE:
!        Calculate the reciver clock corrrection using the
!      presudorange. It is used to calculate the precise satellite coordinate
!
! INPUTS:
!       Obsweek        Observation week
!       Obssec           Observation seconds in a week
!       ObsData         Observation data in a specific station
!       Coor                Coordinate of the station, unit in meter
!       Rotation          Rotation matirx(3*3) from XYZ rectangular coordinate system
!                               to station fixed coordinate system
!      ZHD                  Zenith Hydrostatic Delay, unit in meter
!      ZWD                 Zenith Wet Delay, unit in meter
!      Lat                    Latitude, unit in degree
!      Hgt                   Longitude, unit in degree
!      DOY                  Day of year
!
! OUTPUT:
!        Rec_Clk         receiver clock correction, unit in second
!  
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! =================End of header===============

subroutine Cal_Rec_Clk2(k, Obsweek,Obssec,ObsData,Coor, Rotation,  ZHD,ZWD,Lat, Hgt, DOY, PRNOUT, PRNOUTn, Rec_Clk)
use MOD_constant
use MOD_ObsData
use MOD_Glo_Fre
use MOD_Var
use MOD_FileID
use MOD_ESC
use MOD_CycleSlip
use MOD_STA
use MOD_NavData
implicit none
    ! Intent in
    type(type_ObsData) :: ObsData
    integer Obsweek
    real(8) :: Obssec
    real(8) :: Coor(3)
    real(8) :: Rotation(3,3)
    real(8) :: ZHD, ZWD, lat, Hgt
    integer :: DOY
    ! Intent out
    real(8) :: Rec_Clk, Rec_Clk_PRN(MaxPRN)
    real(8) :: P(MaxPRN)
    
    integer :: i, j,jj, PRN
    real(8) :: KK,freq
    character(1) :: System
    real(8) :: P1, P2, C1, C2, Range, L1, L2, Phase, Range_prn(MaxPRN), Phase_prn(MaxPRN), dOMC(MaxPRN), medOMC
    integer :: RefSat=0, k
    real(8) :: Sat_Coor0(3), Sat_Coor(3), Sat_Vel(3), Sat_XYZ(3),s, t1
    real(8) :: Rela, Ele, maxEle, EleELe(MaxPRN)
    real(8) :: Sat_Clk
    real(8) :: STD, Mfh, Mfw
    integer :: n
    real(8) :: maxV,toe, LimClk
    integer :: maxL, PRNPRN(MaxPRN), PRNOUT(10),PRNOUTn
    
    Rec_Clk=0.0D0
    Rec_Clk_prn=0.d0
    KK=0.0d0
    n=0
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
        else   !  other system
            cycle
        end if
        if (SatSelected(PRN)==0) cycle

        if ((System=="G") .or. (System=='J')) then   ! GPS or QZSS
            f1=10.23d6*154d0
            f2=10.23d6*120d0
        elseif (System=="R") then   ! GLONASS
            freq=Fre_Chann(PRN-GNum)
            f1=(1602.0d0+freq*0.5625d0)*1.0D6   ! f1=(1602.0d0+K*0.5625d0)*1.0d6
            f2=(1246.0d0+freq*0.4375d0)*1.0D6
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
            L1=ObsData%L1(i)
            L2=ObsData%L2(i)
        elseif (freq_comb=='L1L3') then
            P1=ObsData%P1(i)
            P2=ObsData%P3(i)
            if ((P1==0.d0) .and. (C1/=0.d0)) P1=C1+DCB(PRN)*c
            if ((P1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=P1- NavData(PRN)%Nav(2)%TGD(1)*c   ! P1(Y)
            if ((P1==0.d0) .and. (C1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=C1- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL1CA*c ! L1C/A
            if ((P2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=P2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL5Q5*c  ! L5Q5
            L1=ObsData%L1(i)
            L2=ObsData%L3(i)
        elseif (freq_comb=='L2L3') then
            P1=ObsData%P2(i)
            P2=ObsData%P3(i)
            if ((P1==0.d0) .and. (C2/=0.d0)) P1=C2
            if ((P1/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=P1- (f1/f2)**2*NavData(PRN)%Nav(2)%TGD(1)*c  ! P2(Y)
            if ((P1==0.d0) .and. (C2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P1=C2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL2C*c   ! L2C
            if ((P2/=0.d0) .and. (IsCNAV) .and. (System=="G")) P2=P2- NavData(PRN)%Nav(2)%TGD(1)*c+NavData(PRN)%Nav(2)%ISCL5Q5*c  ! L5Q5
            L1=ObsData%L2(i)
            L2=ObsData%L3(i)
        end if
        if ((P1 /=0.0d0) .and. (P2 /=0.0d0)) then
            Range=(f1*f1*P1-f2*f2*P2)/(f1+f2)/(f1-f2)
        else
            Range=0.d0
        end if   ! if ((P1 /=0.0) .and. (P2 /=0)) then

        if ((L1 /=0.0d0) .and. (L2 /=0.0d0)) then
            ! Ionospheric-free combination, in meter
            Phase=(f1*L1-f2*L2)/(f1+f2)/(f1-f2)*c    ! IF=(f1*L1-f2*L2)/(f1^2-f2^2)*c
        else
            Phase=0.d0
        end if
        if (index(ObsCombine,"PC")==0) then
            if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
                Range=P1
                Phase=L1*c/f1
            elseif ( ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) .and. (freq_comb=='L2L3')) then
                Range=P1
                Phase=L1*c/f1
            else
                Range=P2
                Phase=L2*c/f2
            end if
        end if
        if (Range==0.d0) cycle
            
        ! Calculate the satellite coordinate and the distance between satellite and reciver
        t1=Range/c
        call Cal_Sat_PosClk(System, Obsweek,Obssec+ObsData%Clk_Bias, PRN, Coor, t1, .true., Sat_Coor0, Sat_Coor, Sat_Vel, Sat_Clk, s, Rela,toe)
        if ( all(abs(Sat_Coor-9999.0d0)<0.1d0) ) cycle   ! If no data of this PRN
        if (dabs(Sat_Clk-9999.0d0)<0.1d0) cycle  ! If no this Satellite clock data
        ! Add the Equivalent Satllite Clock file
        if (dabs(ESC(PRN)-9999.d0)<0.1d0 ) cycle
        Sat_Clk=Sat_Clk-ESC(PRN)/c  ! 加钟差改正

        Sat_XYZ=MATMUL(Rotation,(Sat_Coor-Coor))
        Ele=dasind(Sat_XYZ(3)/dsqrt(DOT_PRODUCT(Sat_XYZ,Sat_XYZ)))
        if (Ele>LimEle) then
            !k=1.d0 ! 2*dsind(Ele)
        else
            cycle   ! The cutoff satellite elevation
        end if

        call Map_NMF(Ele, Lat, Hgt, DOY, MFh, MFw)
        STD=ZHD*Mfh+ZWD*Mfw
        KK=KK+1.d0 ! The coefficient
        n=n+1
        P(n)=1.d0
        PRNPRN(n)=PRN
        Range_prn(N)=Range
        Phase_prn(N)=Phase
        EleEle(n)=Ele
        ! The recvier clock corrction for each satellite
        if ( (ObsType=='X71') .or. (ObsType=='x71') ) then
            if (ObsData%STD(i)/=0.d0) then
                STD=0.d0
            end if
            rela=0.d0
        end if
        Rec_Clk_prn(n)=(Range - s -STD+Rela ) /c +Sat_Clk
!            if ( (i>1) .and. (dabs(Rec_Clk - Rec_Clk_prn)>1.0D-7) .and. (dabs(Rec_Clk_prn)>1.0D-6)) then
        Rec_Clk=Rec_Clk + 1.d0/KK*(Rec_Clk_prn(n) - Rec_Clk)
!                write(unit=LogID,fmt='(A7,I3,3F15.3,3F8.2,2F10.3)') 'PRN',PRN, Range, s,Sat_Clk*c, &
!                                        STD, rela, 0.d0, Rec_Clk_prn(n)*c,Rec_Clk*c
                        
!!            else ! If the diff is less than 1D-7 s or 
!                    !the clock correction is less than 1D-6 s, then we can ignore it.
!!                exit
!!            end if  ! if ((abs(Rec_Clk - R
!        else   ! GLONASS or other system
!            cycle
!        end if   ! if (System=="G")
    end do  ! do i=1,ObsData%PRNS
    
    PRNOUTn=0
    if (n<1) return
    PRNOUT=0
    if (k<=STA%FixNum) then
        LimClk=2.d-7  ! If station fixed, we assume the tolerance of outlier should be strict. However, there is ISB
    else
        LimClk=3.d-7
    end if
    do while (.true.)
        maxV=maxval(dabs(Rec_Clk_prn(1:N)-Rec_Clk))
        maxL=maxloc(dabs(Rec_Clk_prn(1:N)-Rec_Clk),dim=1)
        if  (dabs(maxV)>LimClk) then
            PRNOUTn=PRNOUTn+1
            if (PRNOUTn>=10) return
            Rec_Clk=(Rec_Clk*KK-Rec_Clk_prn(maxL)*P(maxL))/(KK-P(maxL))
            KK=KK-P(maxL)
            PRNout(PRNOUTn)=PRNPRN(maxL)
            write(LogID,'(3X,I3,A4,A36)') PRNPRN(maxL), 'PRN',  '的O-C与均值差异大于1.5d-6s，剔除。'
            do j=1,PRNOUTn
                do jj=1,n
                    if (PRNOut(j)==PRNPRN(jj)) then
                        Rec_Clk_prn(jj)=Rec_Clk
                        exit
                    end if
                end do
            end do
            
        else
            exit
        end if
    end do

    maxEle=0.d0
    if (index(CSmethod,"DD")==0) then
        return
    end if
    ! Cycle Slip detect for single frequency
!    if ( (index(ObsCombine,"PC")==0) .and. (Combination(2) .or. (Var_smooth=="y") .or. (Var_smooth=="Y")) ) then
        do i=1,N
            PRN=PRNPRN(i)
            if (Phase_prn(i)/=0.d0) then
                if ((index(Clk,"CLK")/=0) .and. (index(Orbit,"SP3")/=0)) then ! It's only for stations that precise coordinates are known
                    CycleSlip(k)%OMC(PRN,4)=Rec_Clk_prn(i)*c-Range_prn(i)+Phase_prn(i)-Rec_Clk
                else
                    CycleSlip(k)%OMC(PRN,4)=Phase_prn(i)
                end if
            else
                CycleSlip(k)%OMC(PRN,4)=0.d0
            end if
            if ( (EleEle(i)>maxEle) .and. (all(CycleSlip(k)%OMC(PRN,3:4)/=0.d0)) .and. (index(Clk,"CLK")/=0) .and. (index(Orbit,"SP3")/=0)) then
                maxEle=EleEle(i)
                RefSat=PRN   ! Select a reference satellite
            elseif ( (EleEle(i)>maxEle) .and. (all(CycleSlip(k)%OMC(PRN,1:4)/=0.d0)) ) then
                maxEle=EleEle(i)
                RefSat=PRN   ! Select a reference satellite
            end if
        end do

        dOMC=99.d0
        do i=1,N  ! Satellite difference OMC
            PRN=PRNPRN(i)
            if ( (index(Clk,"CLK")/=0) .and. (index(Orbit,"SP3")/=0) ) then
                if (all(CycleSlip(k)%OMC(PRN,3:4)/=0.d0)) then
                    dOMC(i)= CycleSlip(k)%OMC(PRN,4)-CycleSlip(k)%OMC(PRN,3) &
                                        -CycleSlip(k)%OMC(RefSat,4)+CycleSlip(k)%OMC(RefSat,3)
                elseif (CycleSlip(k)%OMC(PRN,2)/=0.d0) then
                    if ((CycleSlip(k)%OMC(RefSat,4)/=0.d0) .and.  (CycleSlip(k)%Slip(RefSat)==0)) then
                    dOMC(i)= CycleSlip(k)%OMC(PRN,4)-CycleSlip(k)%OMC(PRN,2) &
                                      -CycleSlip(k)%OMC(RefSat,4)+CycleSlip(k)%OMC(RefSat,2)
                    end if
                else
                    dOMC(i)=99.d0
                end if
!            elseif (all(CycleSlip(k)%OMC(PRN,2:4)/=0.d0))  then
!                ! Quadratic polyfit
            end if
        end do

        call median(dOMC,n,medOMC)
        if (abs(medOMC)<0.01d0) medOMC=0.d0 ! medOMC=0.d0 if no cycle slip happens
        CycleSlip(k)%CScount=0
        j=N
        do i=1,N
            PRN=PRNPRN(i)
            CycleSlip(k)%Slip(PRN)=0 ! Initial the cycle slip as 0
            if ( (abs(dOMC(i)-medOMC)>0.07d0) .and. (abs(dOMC(i)-99.d0)>=1.d-8) ) then
                write(CSID,'(2I3,F6.2,F10.2,I3)') k, PRN, EleEle(i), dOMC(i)-medOMC, 1
                CycleSlip(k)%Slip(PRN)=1
                CycleSlip(k)%CScount=CycleSlip(k)%CScount+1
            elseif (abs(dOMC(i)-99.d0)<1.d-8) then
                write(CSID,'(2I3,F6.2,F10.2,I3)') k, PRN, EleEle(i), dOMC(i)-medOMC, 2
                CycleSlip(k)%Slip(PRN)=1   ! new satellite appear
                j=j-1
            end if
        end do
        if ((CycleSlip(k)%CScount*2>=j) .or. (j*2<N)) then ! If half of the satellites cycle slips, this method is no longer valid
            CycleSlip(k)%CScount=99  ! 卫星太少或周跳太多时设置99没错，只能怪后面GF探不出来
        end if

        ! shift OMC
        CycleSlip(k)%OMC=eoshift(CycleSlip(k)%OMC,shift=1,dim=2,boundary=0.d0)

!    end if

    return
end subroutine

subroutine median(a,N,med)
implicit none
integer :: N
real(8) :: a(N), med
real(8) :: b(N), x
integer :: i, j, k, NN
   
   NN=N
   med=0.d0
   do i=1,N
       if (a(i)==99.d0) NN=NN-1
   end do

   do i=1,N
       k=0
       if (a(i)==99.d0) cycle
!       call random_number(x)
!       a(i)=a(i)+x/1.e9
       do j=1,N
           if ((a(i)>a(j)) .and. (a(j)/=99.d0)) then
               k=k+1
            end if
        end do
        if (k==int(NN/2)) then
            med=a(i)
            exit
        end if
   end do
   return
end subroutine
    
    

