!   =========Preprocess=======
!
!   Data preprocess , including cycle slip detection and code smooth
!         1. When doing the cycle slip detection, we use the LG combination and MW combination.
!         2. When doing the code smooth, the data is devided into different arcs according to the 
!             cycle slip (step1). When a cycle slip happens in a  PRN, an arc is formed, and the approximate
!             ambiguity can be calculated in this arc (eg. the mean difference between all the accepted code 
!            and phase measurements in this arc.See Bernese GPS Software Version 5.0, P94).
!
! Written by Yize Zhang
! ========End of header=======

subroutine PreProcess
use MOD_constant  ! Speed of light
use MOD_FileID
use MOD_SP3Data
use MOD_ClkData
use MOD_ObsHead
use MOD_GLO_Fre
use MOD_STA
use MOD_ObsData
use MOD_CycleSlip
use MOD_Var
implicit none
    type(type_ObsData) ObsData
    integer :: epoch, PRN, i, k, Kk
    character(1) :: System
     integer :: error
     logical :: flag=.true.
    character(2) :: str_PRN
    integer :: ObsWeek
    real(8) :: ObsSec
    real(8) :: C1, P1, P2, L1, L2, Range, Phase
    integer(1) :: LLI1, LLI2
    integer :: AmbEpo(2,91), n
    real(8) ::  Amb(91),tempamb

    real(8) :: Lat,Lon,Hgt, Rotation(3,3)
    real(8) :: Sat_Coor0(3), Sat_Coor(3), Sat_Vel(3), AppCoor(3) , s, t1, Rela, Sat_Clk,SAT_XYZ(3), ELe
    
    character(len=81) :: line
    
   do i=1,91
        AmbID(i)=FileID_Mark
        FileID_Mark=FileID_Mark+1
        PRN=i
        write(str_PRN,"(I2)") PRN
        open(unit=AmbID(i), file="amb_"//str_PRN//".txt",action="readwrite")
    end do
    if (STA%Num>1) then
        write(*,*) "sation number greater than 1, can't using SMT method in this program."
        pause 
        stop
    end if
    
    ! Start data preprocess
    epoch=0
!    AppCoor=ObsHead%AppCoor
!    call XYZ2BLH(AppCoor(1),AppCoor(2),AppCoor(3),Lat,Lon,Hgt)  ! B,L,H
!    ! Rotation matrix from XYZ rectngular coordinate system to station fixed coordinate system
!    Rotation=reshape((/ -dsind(Lat)*dcosd(Lon),  -dsind(Lon),  dcosd(Lat)*dcosd(Lon),  &
!          -dsind(Lat)*dsind(Lon), dcosd(Lon) ,dcosd(Lat)*dsind(Lon),  dcosd(Lat) ,0.d0, dsind(Lat)/), (/3,3/))
!    
    do while(flag) 
        epoch=epoch+1
        write(*,*) "PreProcess:",epoch
        if (epoch==1) then
            Obsweek=ObsHead(1)%GPSweek
            Obssec=ObsHead(1)%GPSsec
            if (GPSweek_st/=0) then
                Obsweek=GPSweek_st
                Obssec=GPSsec_st
            end if
        else
            Obssec=Obssec+Interval
            if (Obssec>=604800.0d0) then
                Obsweek=Obsweek+1
                Obssec=Obssec - 604800.d0
            end if
        end if
        
        if ( (Obsweek-GPSweek_st)*604800.0d0+Obssec-GPSsec_st<0.d0 ) cycle
        if ( (Obsweek-GPSweek_end)*604800.0d0+Obssec-GPSsec_end>0.d0 ) exit
        

!        if (index(CLK,"CLK")/=0) then   ! If Precise Clock
!            ! Detect whether to read 5-min interval clock data
!            do while ( ( (Obsweek-ClkData%GPSweek(2))*604800.0d0+Obssec-ClkData%GPSsec(2) ) > 0.0d0)
!                call ReadClkData()
!                read(unit=ClkID,fmt="(2X)",iostat=error)  ! Check that if reach the end of the clock file
!                backspace(ClkID)
!                if (error /=0) exit
!            end do
!            if ( (( (Obsweek-ClkData%GPSweek(2))*604800.0d0+Obssec-ClkData%GPSsec(2) ) > 300.0d0)  .or. &
!                  (( (Obsweek-ClkData%GPSweek(1))*604800.0d0+Obssec-ClkData%GPSsec(1) ) < -300.0d0) ) then
!                  write(*,*) "Observation time exceeds the Clk data time 5 minutes, please check!"
!            end if
!        end if
        if (index(Orbit,"SP3")/=0) then   ! If Precise Orbit
             ! Detect whether to read 15-min interval SP3 data
            do while ( ( (Obsweek-SP3Data%GPSweek(6))*604800.0d0+Obssec-SP3Data%GPSsec(6) ) > 0.0d0)
                call ReadSP3Data(1)
                read(unit=SP3ID,fmt="(2X)",iostat=error)  ! Check that if reach the end of the SP3 file
                if (error /=0) exit
                backspace(SP3ID)
            end do
            if ( (( (Obsweek-SP3Data%GPSweek(10))*604800.0d0+Obssec-SP3Data%GPSsec(10) ) > 40.0d0)  .or. &
                  (( (Obsweek-SP3Data%GPSweek(1))*604800.0d0+Obssec-SP3Data%GPSsec(1) ) < -40.0d0) ) then
                  write(*,*) "Observation time exceeds the SP3 data time 15mins, please check!"
            end if
        end if

        ! *******For each station********
        do k=1,STA%Num
            AppCoor=STA%STA(k)%Coor(1:3)
            Rotation=STA%STA(k)%Rotation

            call ReadObsData(Obsweek,Obssec,ObsData,k)
!            CycleSlip(k)%dT=CycleSlip(k)%dT+Interval
            read(unit=ObsID(k),fmt="(2X)",iostat=error)
            backspace(ObsID(k))
            if (error /=0) flag=.false.  ! Check that if reach the end of observation file.
            if (ObsData%Flag /=0) cycle
            
        
            ! Cycle slip detect
            do i=1,ObsData%PRNS
                PRN=ObsData%PRN(i)
                System=ObsData%System(i)
                if (System=="G") then  !! GPS
                    if (.not. SystemUsed(1)) cycle
                else if (System=="R") then  ! GLONASS
                    PRN=PRN+32
                    if (.not. SystemUsed(2)) cycle
                else if (System=="C") then  ! COMPASS
                    PRN=PRN+56
                    if (.not. SystemUsed(3)) cycle
                else   ! Galileo or other system
                    cycle
                end if
            
                ! Range
                C1=ObsData%C1(i)
                if (freq_comb=='L1L2') then
                    P1=ObsData%P1(i)
                    P2=ObsData%P2(i)
                    L1=ObsData%L1(i)
                    L2=ObsData%L2(i)
                    LLI1=ObsData%LLI1(i)
                    LLI2=ObsData%LLI2(i)
                elseif (freq_comb=='L1L3') then
                    P1=ObsData%P1(i)
                    P2=ObsData%P3(i)
                    L1=ObsData%L1(i)
                    L2=ObsData%L3(i)
                    LLI1=ObsData%LLI1(i)
                    LLI2=ObsData%LLI3(i)
                elseif (freq_comb=='L2L3') then
                    P1=ObsData%P2(i)
                    P2=ObsData%P3(i)
                    L1=ObsData%L2(i)
                    L2=ObsData%L3(i)
                    LLI1=ObsData%LLI2(i)
                    LLI2=ObsData%LLI3(i)
                end if
                if (System=="G") then   ! GPS
                    f1=10.23d6*154d0
                    f2=10.23d6*120d0
                elseif (System=="R") then   ! GLONASS
                    Kk=Fre_Chann(PRN-32)
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
                end if
                if ((P1 /=0.0d0) .and. (P2 /=0.0d0)) then
                    Range=(f1*f1*P1-f2*f2*P2)/(f1+f2)/(f1-f2) ! Ionospheric-free combination
                    t1=Range/c
                else if ((C1 /=0.0d0) .and. (P2 /=0.0d0)) then
                    P1=C1
                    Range=(f1*f1*C1-f2*f2*P2)/(f1+f2)/(f1-f2)
                    t1=Range/c
                else
                    cycle
                end if   ! if ((P1 /=0.0) .and. (P2 /=0))

                ! ***********Cycle Slip preparation ***********
                Sat_Coor=9999.d0
                if (index(Orbit,"SP3")/=0) then   ! If Precise Orbit
                    if ((System=="G") .or. (System=="R"))then   ! GPS or GLONASS
                        call Cal_Sat_Pos_sp3(Obsweek, Obssec+ObsData%Clk_Bias, PRN, AppCoor, t1, .true. , Sat_Coor0,  Sat_Coor, Sat_Vel, s, Rela)
                    elseif (System=="C") then   ! COMPASS
                        call Cal_Sat_Pos_sp3(Obsweek, Obssec+ObsData%Clk_Bias-14.d0, PRN-56, AppCoor, t1, .true. , Sat_Coor0,  Sat_Coor, Sat_Vel, s, Rela)
                    end if
                elseif  (index(Orbit,"BRD")/=0) then
                    if (System=="G") then   ! GPS
                        call Cal_Sat_Pos_n(System, Obsweek, Obssec+ObsData%Clk_Bias, PRN, AppCoor, t1 , Sat_Coor, Sat_Vel, s, Rela)
                    elseif (System=="R") then   ! GLONASS
                        call Cal_Sat_Pos_g(Obsweek, Obssec+ObsData%Clk_Bias-leap_sec, PRN-32, AppCoor, t1,  Sat_Coor, Sat_Vel, s, Rela, Sat_Clk)
                    elseif (System=="C") then   ! COMPASS
                        call Cal_Sat_Pos_c(Obsweek, Obssec+ObsData%Clk_Bias-14.d0, PRN-56, AppCoor, t1, Sat_Coor, Sat_Vel, s, Rela, Sat_Clk)
                    end if
                end if
                if ( all(dabs(Sat_Coor-9999.0d0)<0.1d0) ) cycle   ! If no data of this PRN
!                if (dabs(Sat_Clk-9999.d0)<0.1d0 ) cycle
            
                Sat_XYZ=MATMUL(Rotation,(Sat_Coor-AppCoor))    ! From XYZ to NEU
                Ele=dasind(Sat_XYZ(3)/dsqrt(DOT_PRODUCT(Sat_XYZ,Sat_XYZ)))   ! Satellite elevation
                ! ---------------- End of Cycle Slip preparation ----------------

                call Cycle_Slip_Detect(k, Obsweek, Obssec, P1, P2, L1, L2, LLI1, LLI2, PRN, Ele, CycleSlip(k)%Slip(PRN))

                if (CycleSlip(k)%Slip(PRN)==0) then  ! no cycle slip happens
                    if ( (L1 /=0.d0) .and. (L2 /=0.d0)) then
                        Phase=(f1*L1-f2*L2)/(f1+f2)/(f1-f2)*c
                        tempamb=Range-Phase   ! in distance
                        if (Amb(PRN)==0.d0) then    ! ambiguity in L1
                            AmbEpo(1,PRN)=epoch    ! Start of the arc
                            AmbEpo(2,PRN)=epoch   ! End of the arc
                        else !if  ( dabs(tempamb(1)-Amb1(PRN))<0.01d0) then
                            AmbEpo(2,PRN)=epoch   ! End of the arc
    !                    else
    !                         write(AmbID(1,PRN),'(I4,2I6,F19.2,A20)') PRN, Amb1Epo(1,PRN), Amb1Epo(2,PRN), Amb1(PRN),'L1 diff>20 cycle'
    !                         Amb1Epo(1,PRN)=epoch   ! Start of the arc   
    !                         Amb1Epo(2,PRN)=epoch   ! End of the arc
                        end if
                        n=AmbEpo(2,PRN)-AmbEpo(1,PRN)+1
                        Amb(PRN)=(1.d0-1.d0/n)*Amb(PRN)+1.d0/n*tempamb  ! in distance
                    end if
                elseif (CycleSlip(k)%Slip(PRN)==1) then  ! If a cycle slip in L1 or L2
                    if (Amb(PRN)/=0.d0) then
                        write(AmbID(PRN),'(I4,2I6,F19.3,A20)') PRN, AmbEpo(1,PRN), AmbEpo(2,PRN), Amb(PRN), 'cycle slip'     
                    end if
                    AmbEpo(1,PRN)=epoch   ! Start of the arc
                    AmbEpo(2,PRN)=epoch   ! End of the arc
                    if ( (L1 /=0.d0) .and. (L2 /=0.d0)) then
                        Phase=(f1*L1-f2*L2)/(f1+f2)/(f1-f2)*c 
                        Amb(PRN)=Range-Phase   ! in distance
                    else
                        Amb(PRN)=0.d0
                    end if
                end if    !  if (Slip==0) then 
            end do    ! do i=1,ObsData%PRNS

        end do  ! do k=1,STA%Num
        ! *******End of for each station********
    end do   ! do while(flag)

    do i=1,91
        PRN=i
         write(AmbID(PRN),'(I4,2I6,F19.3,A20)') PRN, AmbEpo(1,PRN), AmbEpo(2,PRN), Amb(PRN), 'end of '
         rewind(AmbID(PRN))
    end do

    do i=1,STA%Num
        rewind(ObsID(i))
        do while(.true.)    ! Read the header of the obsfile for each station
            read(unit=ObsID(i),fmt="(A)") line
            if (index(line(61:80),"END OF HEADER") /= 0) then
                exit
            end if
        end do
    end do
    if (index(Orbit,"SP3")/=0) then   ! If Precise Orbit
        rewind(SP3ID)
        i=1
        do while(.true.)
            read(unit=SP3ID,fmt="(A80)",iostat=error) line
            i=i+1
            if (i==23) exit
            if (error /=0) exit
        end do
        call ReadSP3Data(10)
    end if
!    if (index(Clk,"CLK")/=0) then   ! If Precise Clock
!        rewind(ClkID)
!        do while(.true.)    ! Read the header of the clkfile
!            read(unit=ClkID,fmt="(A)") line
!            if (index(line(61:80),"END OF HEADER") /= 0) then
!                exit
!            end if
!        end do
!        call ReadClkData
!        call ReadClkData
!    end if  !  if (index(Clk,"CLK")/=0) then
    
end subroutine
