! ===========================================
!  This is the main processing  part of SPP and PPP with
! augmentation information.
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ===============  End of Header  ==============

subroutine Process_Corr
use MOD_constant
use MOD_FileID
use MOD_FileDir
use MOD_NavHead
use MOD_NavData
use MOD_ObsHead
use MOD_ObsData
use MOD_GLO_Fre
use MOD_SP3Data
use MOD_ClkData
use MOD_CycleSlip
use MOD_Ant
use MOD_EOP
use MOD_Rotation
use MOD_Sta
use MOD_Res
use MOD_Var
use MOD_ESC
use MOD_NavData
implicit none
    character(200) :: pcorfile2, orbcorrfile2, zonecorrfile2
    real(8) :: Lat, Lon, Hgt
    real(8) :: ZHD, ZWD, Factor
    logical :: flag, Ad_Flag, alive
    integer:: Obsweek, TropWeek, epoch, EpochUsed
    real(8) :: Obssec, TropSow
    integer(2) :: error
    type(type_ObsData) ObsData
    ! CRS to TRS
    real(8) :: MJD,kmjd, Xp,Yp,dUT1,DX00,DY00
    real(8) :: TT
    real(8) :: FHR, dx_SolTide(3), dx_Ocean(3), dx_CMC(3), dx_pole(3)
    integer :: i, N, NN, Num,j, k, freq, GLOIFB=0
    integer :: PRN, PRN_S, PRNPRN(MaxPRN*2), PRNOUT(10), PRNOUTn=0, ObsNum
    character(1) :: System
    real(8) :: P1,P2,P3, C1, C2, L1,L2,Range, Phase
    integer(1) :: LLI1, LLI2
    real(8) ::Amb(SatNum,STA%Num)
    real(8) :: RecPCO(3),SatPCO(3)
    real(8) :: Azi, Ele_Sat, SatPCV(3)=0.d0, StaPCV(3)=0.d0
    real(8) :: AppCoor(3),AppCoor1(3),AppCoor2(3),Coor(3), BLH(3), Sat_Coor0(3),  Sat_Coor(3), Sat_Vel(3), s, t1, Rela, Sat_XYZ(3)
    real(8) :: STD, Ion=0.d0, Ion1=0.d0
    real(8) :: Rec_Clk, Sat_Clk, Ele, EleEle(MaxPRN*2)
    integer :: CuParaNum
    real(8) ::  Windup_previous(SatNum,STA%Num), dx_windup,temp
    real(8) :: P, PP(MaxPRN*2), A(MaxPRN*2,ParaNum+SatNum),  B(MaxPRN,4),L(MaxPRN*2), Nbb(ParaNum+SatNum,ParaNum+SatNum)
    real(8) ::  U(ParaNum+SatNum), invN(ParaNum+SatNum,ParaNum+SatNum), invN2(4,4), dx(ParaNum+SatNum)
    character(2) :: Code(MaxPRN*2)
    real(8) :: maxV, outlier
    integer :: maxL
    real(8) :: V(40), sigma0, PDOP=0.d0, HDOP, VDOP, ISB(2)
    real(8) :: NEU(3),  Mean_Coor(3), RMS(3), Mean_NEU(3)
    real(8) :: UERE=0.d0,toe
    character(100) :: line

    integer(kind=1)  :: I1outerr,t
    real(8) :: iono_14para(14)    
    real(8) ::Phi, Lam, RcIono,RElevation, RAzimuth,VDelay ,IGP_Delay0(320) , IGP_Delay_SiGma0(320)
    integer(4) :: ierr
    
    
    ! *******For each station********
    do k=1,STA%Num

    flag=.true.
    epoch=0
    EpochUsed=0
    Rec_Clk=0.d0
    Windup_previous=0.d0
    Mean_Coor=0.d0
    Mean_NEU=0.d0
    Amb=0.d0
    Nbb=0.d0
    InvN=0.d0
    U=0.d0
    dx=0.d0

    ! %%%%%% open augmentation files %%%%%%%%%%%% 
    ! Open LC correction file
!    open(101,file=trim(OutDir)//"Resi_O_C_"//STA%STA(1)%Name//'.txt', action='read')
    if (proc_mod==3) then
        zonecorrfile=trim(zonecorrfile)//str_ymd
        inquire(file=zonecorrfile,exist=alive)
        if (.not. alive) then
            write(*,*) "zonecorr file: """//trim(zonecorrfile)//""" doesn't exist!"
            pause
            stop
        end if
        ResO_CID=FileID_Mark
        FileID_Mark=FileID_Mark+1 
        open(ResO_CID,file=zonecorrfile,action='read',err=500)        
        Res(1)%sow=-1000.d0
    end if

     ! open the Equivalent Satllite Clock file
!    open(PcorID,file='C:\Users\yizezhang\Desktop\PCLC\ESC\2012183\Coor_2012183_PCHatch.txt', action='read')
!    do while(.true.)
!        read(PcorID,'(A)', iostat=error) line
!        if (index(line,"++Start of result") /= 0)  exit
!    end do
    ! Open PCOR file
    pcorfile2=trim(pcorfile)//str_ymd//'.dat'
    inquire(file=pcorfile2,exist=alive)
    if (.not. alive) then
        write(*,*) "pcor file: """//trim(pcorfile2)//""" doesn't exist!"
        pause
        stop
    end if
    PcorID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(PcorID,file=pcorfile2,action='read',err=300)
    PCORSow=-1000.d0
    ! Open orbit correction file
    if (proc_mod>=2) then
        orbcorrfile2=trim(orbcorrfile)//str_ymd//'.dat'
        inquire(file=orbcorrfile2,exist=alive)
        if (.not. alive) then
            write(*,*) "orbit correction file: """//trim(orbcorrfile2)//""" doesn't exist!"
            pause
            stop
        end if
        OrbCorrID=FileID_Mark
        FileID_Mark=FileID_Mark+1
        open(OrbCorrID,file=OrbCorrfile2,action='read',err=400)
        OrbCorrSow=-1000.d0
    endif
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    do while(flag)
        epoch=epoch+1
        ! Obsservation Time in GPST
        if (epoch==1) then
            Obsweek=ObsHead(1)%GPSweek
            Obssec=ObsHead(1)%GPSsec
            if (GPSweek_st/=0) then
                Obsweek=GPSweek_st
                Obssec=GPSsec_st
            end if
            if ((STA%STA(k)%SKD=="F") .or. (STA%STA(k)%SKD=="f")) then
                Obssec=Obssec-delay0  ! 算OMC时时间提前
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

        if ((STA%STA(k)%SKD=="F") .or. (STA%STA(k)%SKD=="f")) then
            delay=0.d0   ! 算OMC时和实时一样
        else
            delay=delay0
        end if

        ! Calculate the transfor matrix between TRS and CRS
        MJD=Obsweek*7.d0+44244.0d0+Obssec/86400.d0-Leap_sec/86400.d0 ! In UTC
        kmjd=mod(Obssec,86400.d0)/86400.d0
        EOP%Xp=EOP%X(1)+kmjd*(EOP%X(2)-EOP%X(1))
        EOP%Yp=EOP%Y(1)+kmjd*(EOP%Y(2)-EOP%Y(1))
        dUT1=EOP%dUT1(1)+kmjd*(EOP%dUT1(2)-EOP%dUT1(1))
        DX00=EOP%dX(1)+kmjd*(EOP%dX(2)-EOP%dX(1))
        DY00=EOP%dY(1)+kmjd*(EOP%dY(2)-EOP%dY(1))
        call CRS2TRS(Leap_sec, MJD,EOP%Xp,EOP%Yp,dUT1,DX00,DY00,Rota_C2T,Rota_T2C)
        
        ! Calculate the position of Sun and Moon in CRS
        TT=MJD+2400000.5d0+(19.d0+32.184d0)/86400.d0  ! 地球动力学时
        call PLEPH (TT, 11, 3, SunCoor)
        SunCoor=SunCoor*1000.d0   ! Unit in meter
        call PLEPH (TT, 10, 3, MoonCoor)
        MoonCoor=MoonCoor*1000.d0   ! Unit in meter
        
        ! Get Equivalent Satllite Clock
        !   call Get_ESC(Obsweek,ObsSec)
        ! Get PCOR
        if ( (proc_mod==1) .or. (proc_mod==2) ) then
            if (ObsSec-PCORSow+ObsTime-14.d0>=18.d0) then   ! every 18 seconds, 增强导航，每18s取一次
                call Get_PCOR(Obsweek,int((ObsSec+ObsTime-14.d0+0.1d0)/18.d0)*18.d0,IorQ)  ! 3: Q支路, 0 I支路
            end if
        end if
        if (proc_mod==3) then   ! 分区改正验证，每6min取一次
            if (delay==0) then
                if (ObsSec-PCORSow+ObsTime-14.d0>=360.d0) then
                    call Get_PCOR(Obsweek,int((ObsSec+ObsTime-14.d0+0.1d0)/360.d0)*360.d0,IorQ)  ! 实时OMC
                end if
            else
!                if (ObsSec-PCORSow+ObsTime-14.d0>=360.d0+delay) then   ! 算改正数时delay=0, 和实时OMC一样
!                    call Get_PCOR(Obsweek,int((ObsSec+ObsTime-14.d0-delay+0.1d0)/360.d0)*360.d0,IorQ)  ! 90s或3min延迟OMC，取整6min
!                end if
                if (ObsSec-PCORSow+ObsTime-14.d0>=90.d0+delay) then   ! 算改正数时delay=0, 和实时OMC一样
                    call Get_PCOR(Obsweek,int((ObsSec+ObsTime-14.d0+0.1d0)/delay)*delay-90.d0,IorQ)  ! 90s或3min延迟OMC，取往前90s
                end if
!                if (ObsSec-PCORSow+ObsTime-14.d0>=36.d0+delay) then   ! 算改正数时delay=0, 和实时OMC一样
!                    call Get_PCOR(Obsweek,int((ObsSec+ObsTime-14.d0+0.1d0)/delay)*delay-36.d0,IorQ)  ! 90s或3min延迟OMC，取往前36s
!                end if
            end if
        end if
        ! get X37a
        if (proc_mod==2) then
            if (ObsSec-OrbCorrSow+ObsTime-14.d0>=360.d0) then   ! every 6 minutes, 增强导航，每6min取一次
                call Get_OrbCorr(Obsweek,int((ObsSec+ObsTime-14.d0+0.1d0)/360.d0)*360.d0-10.d0,IorQ)  ! 3: Q支路, 0 I支路
            end if
        end if
        if (proc_mod==3) then   ! 分区改正验证，每6min取一次
            if (delay==0) then
                if (ObsSec-OrbCorrSow+ObsTime-14.d0>=360.d0) then
                    call Get_OrbCorr(Obsweek,int((ObsSec+ObsTime-14.d0+0.1d0)/360.d0)*360.d0,IorQ)  !  ! 实时OMC
                end if
            else
!                if (ObsSec-OrbCorrSow+ObsTime-14.d0>=360.d0+delay) then   ! 定位时,算改正数时delay=0, 和实时OMC一样
!                    call Get_OrbCorr(Obsweek,int((ObsSec+ObsTime-14.d0-delay+0.1d0)/360.d0)*360.d0,IorQ)  ! 90s或3min延迟OMC
!                end if
                if (ObsSec-OrbCorrSow+ObsTime-14.d0>=360.d0+delay+10.d0) then   ! 定位时,算改正数时delay=0, 和实时OMC一样
                    call Get_OrbCorr(Obsweek,int((ObsSec+ObsTime-14.d0+0.1d0)/360.d0)*360.d0-10.d0,IorQ)  ! 90s或3min延迟OMC,X37a提前了10s
                end if
            end if
        end if

        ! *******For each station********
!        do k=1,STA%Num
            if ((proc_mod>=1) .and. (proc_mod<=2) .and. (clktype==2)) then
                    if ((mod(ObsSec,3600.d0)+ObsTime-14.d0>1980) .and. (mod(ObsSec,3600.d0)+ObsTime-14.d0<2100)) then
                        cycle  ! 实时，切换星历时不用
                    end if
            elseif ((proc_mod>=1) .and. (proc_mod<=2) .and. (clktype==1)) then
                    if ((mod(ObsSec,3600.d0)+ObsTime-14.d0>0) .and. (mod(ObsSec,3600.d0)+ObsTime-14.d0<30)) then
                        cycle  ! 实时，切换星历时不用
                    end if
            end if

            ! 获取伪距相位改正数 ! get LC corrections
            if (proc_mod==3) then
                if ((STA%STA(k)%SKD=="K") .or. (STA%STA(k)%SKD=="k")) then
                    if (ObsSec-Res(1)%sow+ObsTime-14.d0>=delay) then
!                        if (delay==0.d0) then
!                            call Get_Corr(Obsweek, ObsSec+ObsTime-14.d0-delay)  ! delay=0.d0
!                        else
!                            call Get_Corr(Obsweek,int((ObsSec+ObsTime-14.d0+0.1d0)/delay)*delay)
!                        end if
                        if (delay==0.d0) then
                            if ((mod(ObsSec,3600.d0)+ObsTime-14.d0>1980) .and. (mod(ObsSec,3600.d0)+ObsTime-14.d0<2100)) then
                                cycle  ! 实时，切换星历时不用
                            end if
                            call Get_X38(Obsweek,ObsSec+ObsTime-14.d0)
                        else
!                            if ((mod(ObsSec,3600.d0)+ObsTime-14.d0>36) .and. (mod(ObsSec,3600.d0)+ObsTime-14.d0<108)) then
!                                cycle  ! 实时，切换星历时不用
!                            end if
                            call Get_X38(Obsweek,int((ObsSec+ObsTime-14.d0+0.1d0)/delay)*delay)
                        end if
                    end if
                end if
            end if

            AppCoor=STA%STA(k)%Coor(1:3)+MATMUL( STA%STA(k)%NEU, STA%STA(k)%Rotation )
            ZHD=STA%STA(k)%Trop%ZHD
            ZWD=STA%STA(k)%Trop%ZWD
            Rotation=STA%STA(k)%Rotation
            Lat=STA%STA(k)%BLH(1)
            Lon=STA%STA(k)%BLH(2)
            Hgt=STA%STA(k)%BLH(3)
            
            ! Read obs data at the given GPST
            if ( (obstype=='RNX') .or. (obstype=='rnx') ) then
                call ReadObsData(Obsweek,Obssec,ObsData,k)
            elseif  ( (ObsType=='X71') .or. (ObsType=='x71') ) then
                call ReadX71Data(Obsweek,Obssec,ObsData,k)
            elseif  ( (ObsType=='X11') .or. (ObsType=='x11') ) then
                call ReadX11Data(Obsweek,Obssec,ObsData,k)
            end if
!            CycleSlip(k)%dT=CycleSlip(k)%dT+Interval
            read(unit=ObsID(k),fmt="(2X)",iostat=error)
            backspace(ObsID(k))
            if (error /=0) flag=.false.  ! Check that if reach the end of observation file.
            if (ObsData.Flag /=0) cycle
            
            if (mod(epoch,10)==1) write(*,'(A5,I6,2X,I4,4I3,F5.1)') "Epoch", epoch, ObsData%year,ObsData%mon,&
                           ObsData%day, ObsData%hour, ObsData%min, ObsData%sec

            ! Calculate the appropriate position if unknown
!            if ( any(STA%STA(k)%TrueCoor==0.d0) .or. (Pos_State=='K') ) then
!            if ( any(STA%STA(k)%TrueCoor==0.d0) .or. ((Pos_State=='K') .and. .not.(STA%STA(k)%flag_InitialCoor))  ) then
            if ( any(STA%STA(k)%TrueCoor==0.d0) .or. (Pos_State=='K' .and. any(dabs(STA%STA(k)%TrueCoor-STA%STA(k)%XYZ)>30.d0))  ) then
                call Bancroft(ObsData, k, STA%STA(k)%Coor)
                AppCoor=STA%STA(k)%Coor
                if (any(AppCoor==0.d0)) cycle
                ZHD=STA%STA(k)%Trop%ZHD
                ZWD=STA%STA(k)%Trop%ZWD
                Rotation=STA%STA(k)%Rotation
                Lat=STA%STA(k)%BLH(1)
                Lon=STA%STA(k)%BLH(2)
                Hgt=STA%STA(k)%BLH(3)
            end if
            
            ! 1. Calculate the solid tide correction
            FHR=(ObsData%hour+ObsData%min/60.d0+ObsData%sec/3600.d0)
            dx_SolTide=0.d0
            call Tide(AppCoor,ObsData%year,ObsData%mon,ObsData%day,FHR, &
     &     MATMUL(Rota_C2T,SunCoor(1:3)), MATMUL(Rota_C2T,MoonCoor(1:3)), dx_SolTide)  ! In TRS
            
            if ( (ObsType=='X71') .or. (ObsType=='x71') ) dx_SolTide=0.d0
            AppCoor1=AppCoor +  dx_SolTide    ! For Process_Corr, we don't add ocean load tide and pole tide corrections

             ! ********* Calculate the initial value of reciver clock corrrection *********
             !call Cal_Rec_Clk(Obsweek, Obssec, ObsData,AppCoor1, Rotation, Rec_Clk)
!             Rec_Clk=0.d0
             call Cal_Rec_Clk2(k, Obsweek, Obssec, ObsData,AppCoor1, Rotation, ZHD,ZWD,Lat, Hgt, int_doy, PRNOUT, PRNOUTn, Rec_Clk)
             if ((STA%STA(k)%SKD=="K") .or. (STA%STA(k)%SKD=="k")) then
                 if (ObsData%PRNS - PRNOUTn <=3) then
                     write(unit=LogID,fmt='(A7, I3,A30)')  'Epoch', epoch,'too few satellites, skip.'
                     cycle
                 end if
                 
             else
!                 write(unit=LogID,fmt='(A5,I5,2X, 5I4,F6.1,E15.7)') 'Epoch', epoch, ObsData%year,ObsData%mon,&
!                       ObsData%day, ObsData%hour, ObsData%min, ObsData%sec, Rec_Clk
             end if
             if (ObsHead(k)%Version==2) then
                 write(unit=LogID,fmt='(A5,I5,2X,I4,4I3,F5.1,F15.3,I5,I8,3F10.3)') 'Epoch', epoch, ObsData%year,ObsData%mon,&
                           ObsData%day, ObsData%hour, ObsData%min, ObsData%sec, Rec_Clk*c, Obsweek, int(Obssec), dx_SolTide
            else
                 write(unit=LogID,fmt='(A5,I5,2X,I4,4I3.2,F5.1,F15.3,I5,I8,3F10.3)') 'Epoch', epoch, ObsData%year,ObsData%mon,&
                           ObsData%day, ObsData%hour, ObsData%min, ObsData%sec, Rec_Clk*c, Obsweek, int(Obssec), dx_SolTide
            end if

            ! ========Positioning============
                N=0
                NN=0
                ObsNum=0
                PRNPRN=0
                EleEle=0.d0
                PP=0.d0
                Code=""
                A=0.d0
                L=0.d0
                Range=0.d0
                CuParaNum=0
                write(CSID,"(A6,I5)")  "epoch:", epoch
                if ((STA%STA(k)%SKD=="F") .or. (STA%STA(k)%SKD=="f")) then
                    if (delay0==0.d0) then
                        write(ResO_CID,"(A1,1X,I4,F10.1, 2X, I4,4(I3),F5.1)")  '*', Obsweek,Obssec, ObsData%year, ObsData%mon, ObsData%day, ObsData%hour,ObsData%min, ObsData%sec
                    else
                        write(ResO_CID,"(A1,1X,I4,F10.1, 2X,I4,4(I3),F5.1)")  '*', Obsweek,ceiling(ObsSec/delay0)*delay0, &
                         ObsData%year, ObsData%mon, ObsData%day, ObsData%hour,ObsData%min, ObsData%sec
                    end if
                end if
                do i=1,ObsData%PRNS
                    PRN_S=ObsData%PRN(i)
                    PRN=PRN_S
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
                        100 cycle
                    end if
                    do j=1, PRNOUTn
                        if (PRN==PRNOUT(j)) goto 100
                    end do
                    if (SatSelected(PRN)==0) cycle
                    
                    if ((proc_mod==3) .and. (Res(1)%L((PRN-56)*2)==0.d0)) then
                        cycle
                    end if

                    if ((System=="G") .or. (System=='J')) then   ! GPS/QZSS
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

                    ! Range
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
                        LLI1=ObsData%LLI1(i)
                        LLI2=ObsData%LLI2(i)
                        if (L2==0.d0 .and. ObsData%L2C(i)/=0.d0) then
                            L2=ObsData%L2C(i)
                            LLI2=ObsData%LLI2C(i)
                        end if
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
                        if (L1==0.d0 .and. ObsData%L2C(i)/=0.d0) then
                            L1=ObsData%L2C(i)
                            LLI1=ObsData%LLI2C(i)
                        end if
                    end if
                    if ((P1 /=0.0d0) .and. (P2 /=0.0d0)) then
                        Range=(f1*f1*P1-f2*f2*P2)/(f1+f2)/(f1-f2) ! Ionospheric-free combination
                    else if ((C1 /=0.0d0) .and. (P2 /=0.0d0)) then
                        P1=C1  ! For net difference
                        Range=(f1*f1*C1-f2*f2*P2)/(f1+f2)/(f1-f2)
                    elseif (index(ObsCombine,"PC")==0) then ! If single frequency
                        if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
                            Range=P1
                        elseif ( ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) .and. (freq_comb=='L2L3')) then
                            Range=P1
                        else   ! G2 or (PC without P1 or P2)
                            Range=P2
                        end if
                    else
                        Range=0.d0
                    end if
                    if (Range==0.d0) cycle
                    t1=Range/c
                    
                    ! Receiver PCO
                    if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
                        RecPCO=Ant(GNum+RNum+CNum+NumE+JNum+INum+k)%PCO(1:3,1)
                    elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) then
                        RecPCO=Ant(GNum+RNum+CNum+NumE+JNum+INum+k)%PCO(1:3,2)
                    else
                        RecPCO=(f1*f1*Ant(GNum+RNum+CNum+NumE+JNum+INum+k)%PCO(1:3,1)-f2*f2*Ant(GNum+RNum+CNum+NumE+JNum+INum+k)%PCO(1:3,2))/(f1+f2)/(f1-f2) 
                    end if
                    AppCoor2=AppCoor1 +MATMUL(transpose(Rotation),RecPCO)

                    ! ********Satellite coordinate and clcok ********
                    call Cal_Sat_PosClk(System, Obsweek, Obssec+ObsData%Clk_Bias-Rec_Clk,PRN, AppCoor2, t1, .true., Sat_Coor0,  Sat_Coor, Sat_Vel, Sat_Clk, s, Rela,toe)
                    if ( all(dabs(Sat_Coor-9999.0d0)<0.1d0) ) cycle   ! If no data of this PRN
                    if (dabs(Sat_Clk-9999.d0)<0.1d0 ) cycle
                    ! Add the Equivalent Satllite Clock file
                    if (dabs(ESC(PRN)-9999.d0)<0.1d0 ) cycle
                    if (proc_mod==3) then
                        if (dabs(ESC(PRN)-0.d0)<0.001d0 ) cycle ! 没有钟差改正数时不计算
                    end if
                    Sat_Clk=Sat_Clk-ESC(PRN)/c  ! 加钟差改正
                    if (proc_mod>=2) then  ! 加轨道改正
                        Sat_Coor=Sat_Coor +OrbCorr(PRN,1:3)
                        s=dsqrt(DOT_PRODUCT((Sat_Coor-AppCoor2),(Sat_Coor-AppCoor2)))
                    end if

                    Sat_XYZ=MATMUL(Rotation,(Sat_Coor-AppCoor2))
                    Ele=dasind(Sat_XYZ(3)/dsqrt(DOT_PRODUCT(Sat_XYZ,Sat_XYZ)))
                    Azi=datand(Sat_XYZ(2)/Sat_XYZ(1))    ! Satellite Azimuth
                    if (Sat_XYZ(1)<0.d0) then
                        Azi=Azi+180.d0
                    elseif (Sat_XYZ(2)<0.d0) then
                        Azi=Azi+360.d0
                    end if
                    if (Ele>30.d0) then
                        P=1.d0
                    elseif (Ele>LimEle) then
                        P=2*dsind(Ele)
                    else
                        write(LogID,'(A6,1X,A1,I2,F8.2,A20)')  'PRN', System, PRN_S, Ele, 'elevation too low'
                        cycle
                    end if
                    
                    ! ********Satellite PCO (Phase Center Offset) correction********
                    ! Transformation matrix from Satellite-fixed system to CRS
                    call S2C(SunCoor(1:3), Sat_Coor, Sat_Vel, System, PRN)
                    if ( (ObsType=='RNX') .or. (ObsType=='rnx') ) then
                        if  (Orbit=="SP3") then
                            if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
                                SatPCO=Ant(PRN)%PCO(1:3,1)
                            elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) then
                                SatPCO=Ant(PRN)%PCO(1:3,2)
                            else
                                SatPCO=(f1*f1*Ant(PRN)%PCO(1:3,1)-f2*f2*Ant(PRN)%PCO(1:3,2))/(f1+f2)/(f1-f2)
                            end if
                            !SatPCO(3) = SatPCO(3) - dPCO(PRN)
                            Sat_Coor=Sat_Coor+MATMUL(Rota_C2T, MATMUL(Rota_S2C,SatPCO))  ! Satellite phase center in TRS
                            ! ********Satellite PCV (Phase Center Variation) correction********
                            Ele_Sat=dasind(sind(Ele+90.d0)*6378137.d0/dsqrt(DOT_PRODUCT(Sat_Coor,Sat_Coor)))
                            call PCV_Corr(PRN, Ele_Sat, 0.d0, SatPCV)  ! Satllite PCV, zenith depended
                            if (allocated(Ant(GNum+RNum+CNum+NumE+JNum+INum+k)%PCV)) then
                                call PCV_Corr(GNum+RNum+CNum+NumE+JNum+INum+k, 90.d0-Ele, Azi, StaPCV)  ! Station PCV, zenith and azimuth depended
                            end if
                            if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
                                SatPCV(3)=SatPCV(1)
                                StaPCV(3)=StaPCV(1)
                            elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) then
                                SatPCV(3)=SatPCV(2)    ! For satellite, PCV on L1 and L2 is the same, actually
                                StaPCV(3)=StaPCV(2)
                            else
                                SatPCV(3)=(f1*f1*SatPCV(1)-f2*f2*SatPCV(2))/(f1+f2)/(f1-f2)
                                StaPCV(3)=(f1*f1*StaPCV(1)-f2*f2*StaPCV(2))/(f1+f2)/(f1-f2)
                            end if
                            s=dsqrt(DOT_PRODUCT((Sat_Coor-AppCoor2),(Sat_Coor-AppCoor2))) +SatPCV(3)+StaPCV(3)     ! real distance
                        end if
                    end if
                    
                    ! 相位平滑伪距
!                    if ( Combination(2)==.false. ) then  ! 没有LC组合的时候才用相位平滑伪距
                        if ( (Var_smooth=="y") .or. (Var_smooth=="Y") ) then
                            if (index(Smooth_Method,"Hatch") /=0)  then
                                call Cycle_Slip_Detect(k, Obsweek, Obssec, P1, P2, L1, L2, LLI1, LLI2, PRN, Ele, CycleSlip(k)%Slip(PRN))
                                call Hatch_Filter(PRN, P1, P2, Range, L1*c/f1, L2*c/f2, Ele, k, CycleSlip(k)%Slip(PRN))
                            elseif (index(Smooth_Method,"SMT") /=0) then
                                ! RNX Smoothing
                                call RNXSMT(epoch, PRN, Range, L1*c/f1, L2*c/f2)
                            end if
                        end if
!                    end if

                    ! ************STD(Slant Tropsphere Delay)*************
                    if  (cmap(1:4)=='SAAS') then
                        STA%STA(k)%Trop%map_dry=1.d0/dsind(Ele)
                        STA%STA(k)%Trop%map_dry=STA%STA(k)%Trop%map_wet
                    elseif (cmap(1:3)=='NMF') then
                        call Map_NMF(Ele, Lat, Hgt, int_doy, STA%STA(k)%Trop%map_dry, STA%STA(k)%Trop%map_wet)
                     elseif (cmap(1:3)=='GMF') then
                        call Map_GMF(MJD,Lat*pi/180.d0,Lon*pi/180.d0, Hgt,  (90.d0-Ele)*pi/180.d0, STA%STA(k)%Trop%map_dry, STA%STA(k)%Trop%map_wet)
                    elseif ( (cmap(1:4)=='VMF1') .and. (cdattype(1:4)=='GPT2') )  then
                        call vmf1_ht(STA%STA(k)%Trop%ah,STA%STA(k)%Trop%aw, MJD, Lat*pi/180.d0, Hgt, (90.d0-Ele)*pi/180.d0, STA%STA(k)%Trop%map_dry, STA%STA(k)%Trop%map_wet)
                        ! call vmf1(STA%STA(k)%Trop%ah,STA%STA(k)%Trop%aw, MJD, Lat*pi/180.d0, (90.d0-Ele)*pi/180.d0, STA%STA(k)%Trop%map_dry, STA%STA(k)%Trop%map_wet)
                    else
                        write(*,*) '***ERROR : unknown ZTD mapping  function or wrong combination: ',cmap,cdattype
                        pause
                        stop
                    end if
                    Factor=STA%STA(k)%Trop%map_wet
                    STD=STA%STA(k)%Trop%map_wet*ZWD+STA%STA(k)%Trop%map_dry*ZHD
                    if ( (ObsType=='x71') .or. (ObsType=='X71') ) then
                        if (ObsData%STD(i)/=0.d0) then
                            STD=0.d0
                        end if
                        if (ObsData%rela(i)/=0.d0) then
                            rela=0.d0
                        end if
                    end if
                    if (TropLen==0.d0) then  ! If Trpsphere parameter is not estimated
                        Factor=0.d0
                    end if
                    
                    if (index(ObsCombine,"PC")==0) then   ! If single frequency 
                        if (iontype==1) then
                            if (any(NavHead%Alpha/=0.d0)) then
                                call Klobuchar(Lat, Lon, Ele, Azi, Obssec, NavHead%Alpha,NavHead%Beta,Ion)
                            else
                                Ion=0.d0
                            end if
                        elseif (iontype==2) then
                            call GIM(Lat, Lon, Ele, Azi, ObsWeek, ObsSec,Ion)  ! 基于B1频点
                        elseif (iontype==3) then
                            ! 利用14参计算北斗B1频点的电离层延迟
                            t=int(mod(Obssec,86400.d0)/7200.d0)+1
                            iono_14para=iono_14para_day(t,1:14)
                            call cal_iondelay(int(Obssec),AppCoor2,Sat_Coor,iono_14para, Ion, I1outerr)  ! B1 iono delay
                        elseif (iontype==4) then  !BD 格网
                            ! 利用14参计算北斗B1频点的电离层延迟
                            t=int(mod(Obssec,86400.d0)/7200.d0)+1
                            iono_14para=iono_14para_day(t,1:14)
                            call cal_iondelay(int(Obssec),AppCoor2,Sat_Coor,iono_14para, Ion1, I1outerr)  ! B1 iono delay

                            call Pierce_PHiLam(AppCoor2(1),AppCoor2(2),AppCoor2(3),Sat_Coor(1),Sat_Coor(2),Sat_Coor(3), Phi, Lam, RcIono, &
                                       RElevation, RAzimuth,Ierr)
                            call get_igpdely(int4(Obsweek-1356),Obssec,IGP_Delay0,IGP_Delay_SiGma0,ierr)
                            call IntpoGrid(Phi/PI*180.d0 , Lam/PI*180.d0 , VDelay ,IGP_Delay0 , IGP_Delay_SiGma0,Ierr)
                            Ion = VDelay/RcIono  ! 基于B1频点
                            if (Ion==0.d0) then
                                cycle   ! If no grid data
                                Ion=Ion1
                            end if
                        end if
                    end if
                    if ( (index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0) )then
                        if (iontype<=2) then
                            Ion=Ion*(154.d0*10.23d6)**2/f1**2
                        end if
                        if (P1 /=0.0d0) Range=P1-Ion  ! only use P1 and L1
                        if (L1 /=0.0d0) then
                            if (index(ObsCombine,"G1")/=0) then
                                Phase=(ObsData%L1(i)*c/f1+P1)/2.d0   ! GRAPHIC combination, (P1+L1)/2
                            elseif (index(ObsCombine,"P1")/=0) then
                                Phase=ObsData%L1(i)*c/f1+Ion
                            end if
                        else
                            Phase=0.d0
                        end if
                    elseif ( (index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0) )then
                        ! B2频点的电离层延迟
                        if (iontype<=2) then
                        Ion=Ion*(154.d0*10.23d6)**2/f2**2
                        else
                        Ion=Ion*f1**2/f2**2
                        end if
                        if (P2 /=0.0d0) Range=P2-Ion  ! only use P2 and L2
                        if (L2 /=0.0d0) then
                            if (index(ObsCombine,"G2")/=0) then
                                Phase=(ObsData%L2(i)*c/f2+P2)/2.d0   ! GRAPHIC combination, (P2+L2)/2
                            elseif (index(ObsCombine,"P2")/=0) then
                                Phase=ObsData%L2(i)*c/f2+Ion
                            end if
                        else
                            Phase=0.d0
                        end if
                    elseif ( (index(ObsCombine,"P3")/=0) .or. (index(ObsCombine,"G3")/=0) )then
                        ! B3频点的电离层延迟
                        if (iontype<=2) then
                        Ion=Ion*(154.d0*10.23d6)**2/f2**2
                        else
                        Ion=Ion*f1**2/f2**2
                        end if
                        if (P2 /=0.0d0) Range=P2-Ion  ! only use P3 and L3
                        if (L2/=0.0d0) then   ! L2=ObsData%L3(i), P2=ObsData%P3(i)
                            if (index(ObsCombine,"G3")/=0) then
                                Phase=(ObsData%L3(i)*c/f2+P2)/2.d0   ! GRAPHIC combination, (P3+L3)/2
                            elseif (index(ObsCombine,"P3")/=0) then
                                Phase=ObsData%L3(i)*c/f2+Ion
                            end if
                        else
                            Phase=0.d0
                        end if
                    elseif (index(ObsCombine,"PC")/=0) then
                        if ((L1 /=0.0d0) .and. (L2 /=0.0d0)) then
                            ! Ionospheric-free combination, in meter
                            Phase=(f1*L1-f2*L2)/(f1+f2)/(f1-f2)*c    ! IF=(f1*L1-f2*L2)/(f1^2-f2^2)*c
                        else
                            Phase=0.d0
                        end if
                    end if

                    ! ******************** for range *****************
                    if (Range/=0.d0) then
                        N=N+1
                        NN=NN+1
                        PRNPRN(N)=PRN
                        Code(N)="PC"
                        EleEle(N)=Ele
                        PP(N)=P/sigPC
                        A(N,:)=0.d0
                        A(N,1:4)=(/(AppCoor2-Sat_Coor)/s,  Factor/)*PP(N)
                        A(N,5)=PP(N)
!                        A(N, PRN+ParaNum)=PP(N)
                        B(NN,:)=0.d0
                        B(NN,1:4)=(/(AppCoor2-Sat_Coor)/s, 1.d0/)
                        if ((System=="R") .and. (If_ISB)) then
                            if (IF_IFB) then
                                if (IFB_Mode=='FD') then
                                    A(N,freq+13+INT_SystemUsed(1))=PP(N) ! Frequency depend
                                    GLOIFB=14
                                elseif (IFB_Mode=='LM') then
                                    A(N,5+INT_SystemUsed(1))=PP(N)   ! linear model
                                    A(N,6+INT_SystemUsed(1))=PP(N)*freq
                                    GLOIFB=1
                                end if
                            else
                                A(N,5+INT_SystemUsed(1))=PP(N)
                                GLOIFB=0
                            end if
                        elseif ((System=="C") .and. (If_ISB)) then
                            A(N,5+INT_SystemUsed(1)+GLOIFB+INT_SystemUsed(2))=PP(N)
                        elseif ((System=="E") .and. (If_ISB)) then
                            A(N,5+INT_SystemUsed(1)+GLOIFB+INT_SystemUsed(2)+INT_SystemUsed(3))=PP(N)
                        elseif ((System=="J") .and. (If_ISB)) then
                            A(N,5+INT_SystemUsed(1)+GLOIFB+INT_SystemUsed(2)+INT_SystemUsed(3)+INT_SystemUsed(4))=PP(N)
                        elseif ((System=="I") .and. (If_ISB)) then
                            A(N,5+INT_SystemUsed(1)+GLOIFB+INT_SystemUsed(2)+INT_SystemUsed(3)+INT_SystemUsed(4)+INT_SystemUsed(5))=PP(N)
                        end if
                        if (If_Est_Iono .and. (index(ObsCombine,"P1") /=0 .or. index(ObsCombine,"P2") /=0 .or. index(ObsCombine,"P3") /=0)) then
                            A(N, PRN+ParaNum-SatNum)= PP(N)  ! Ionosphere parameter
                        end if
                        L(N)=(Range-s-STD+Sat_Clk*c-c*Rec_Clk+Rela )
                        if ((STA%STA(k)%SKD=="F") .or. (STA%STA(k)%SKD=="f")) then
                            write(ResO_CID,"(A5,I3,F10.4,I3)")  Code(N), PRNPRN(N), L(N), 0  ! write O-C residuals
                            write(X38ID,'(I4,F10.1,2I3,2X,A4,2X,I3,A4,F10.3,I2)') Obsweek-1356, ObsSec, 5, IorQ, &
                                    STA%STA(k)%NAME , PRN, "PC", L(N), 0
                        end if
                    end if
                     
                    ObsNum=ObsNum+1   
                    ! ================ For Phase =================
                    if ((Combination(2)==.false.) .or. (Phase==0.d0)) then
                        write(unit=LogID,fmt='(A6,1X,A1,I2,F8.2,2F15.3,E15.7,4F8.3,F10.3)') 'PRN',System,PRN_S,Ele, Range, s, Sat_Clk,STD, Ion,rela, ESC(PRN+56), L(N)
                        cycle
                    end if

                    if ( (Var_smooth/="y") .and. (Var_smooth/="Y") ) then
                        call Cycle_Slip_Detect(k, Obsweek, Obssec, P1, P2, L1, L2, LLI1, LLI2, PRN, Ele, CycleSlip(k)%Slip(PRN))
                    end if
                     if ( (arc_epochs/=0) .and. (mod(epoch,arc_epochs)==1) )  then
                        CycleSlip(k)%Slip(PRN)=1
                    end if
                    if (CycleSlip(k)%Slip(PRN)==1) then
                        CuParaNum=CuParaNum+1
                         if (ADmethod=='LS') then
                            call Elimi_Para(Nbb,U,SatNum+ParaNum,PRN+ParaNum)
                        elseif (ADmethod=='KF') then
                            call KF_Change(InvN,dx,SatNum+ParaNum,PRN+ParaNum,'amb')
                        end if
                        Amb(PRN,k)=0.d0
                            if (epoch>1) then
                                write(unit=LogID,fmt='(A6,1X,A1,I2,A24,A5)') 'PRN',System,PRN_S,'cycle slip, at station', STA%STA(k)%Name
                            end if
                    end if
                    
                    ! ***********Phase wind up correction***********
                    call Phase_windup(AppCoor2, Sat_Coor, CycleSlip(k)%Slip(PRN), dx_windup, Windup_previous(PRN,k))
                    dx_windup=0.d0

                    ! *********The initial value of ambiguity*********
                    if (Amb(PRN,k)==0.d0) then
                        Amb(PRN,k)=(Phase-Range)   ! Phase-s-STD+Sat_Clk*c-c*Rec_Clk+Rela -dx_windup*c/(f1+f2) ! In meter   ! randa=c/(f1+f2)
                    end if
                                
                    if (Phase/=0.d0) then
                        N=N+1
                        PRNPRN(N)=PRN
                        Code(N)="LC"
                        EleEle(N)=Ele
                        PP(N)=P/sigLC
                        A(N,:)=0.d0
                        A(N,1:4)=(/(AppCoor2-Sat_Coor)/s,  Factor/)*PP(N)
                        A(N,5)=PP(N)
                        if ((System=="R") .and. (If_ISB)) then
                            if (IF_IFB) then
                                if (IFB_Mode=='FD') then
                                    A(N,freq+13+INT_SystemUsed(1))=PP(N) ! Frequency depend
                                    GLOIFB=14
                                elseif (IFB_Mode=='LM') then
                                    A(N,5+INT_SystemUsed(1))=PP(N)   ! linear model
                                    A(N,6+INT_SystemUsed(1))=PP(N)*freq
                                    GLOIFB=1
                                end if
                            else
                                A(N,5+INT_SystemUsed(1))=PP(N)
                                GLOIFB=0
                            end if
                        elseif ((System=="C") .and. (If_ISB)) then
                            A(N,5+INT_SystemUsed(1)+GLOIFB+INT_SystemUsed(2))=PP(N)
                        elseif ((System=="E") .and. (If_ISB)) then
                            A(N,5+INT_SystemUsed(1)+GLOIFB+INT_SystemUsed(2)+INT_SystemUsed(3))=PP(N)
                        elseif ((System=="J") .and. (If_ISB)) then
                            A(N,5+INT_SystemUsed(1)+GLOIFB+INT_SystemUsed(2)+INT_SystemUsed(3)+INT_SystemUsed(4))=PP(N)
                        elseif ((System=="I") .and. (If_ISB)) then
                            A(N,5+INT_SystemUsed(1)+GLOIFB+INT_SystemUsed(2)+INT_SystemUsed(3)+INT_SystemUsed(4)+INT_SystemUsed(5))=PP(N)
                        end if
                        A(N, PRN+ParaNum)=PP(N)
                        if (If_Est_Iono .and. (index(ObsCombine,"P1") /=0 .or. index(ObsCombine,"P2") /=0 .or. index(ObsCombine,"P3") /=0)) then
                            A(N, PRN+ParaNum-SatNum)= - PP(N)  ! Ionosphere parameter
                        end if
                        L(N)=Phase-s-STD+Sat_Clk*c-c*Rec_Clk+Rela-Amb(PRN,k) -dx_windup*c/(f1+f2) !
                        if ((STA%STA(k)%SKD=="F") .or. (STA%STA(k)%SKD=="f")) then
                            write(ResO_CID,"(A5,I3,F10.4,I3)")  Code(N), PRNPRN(N), L(N), CycleSlip(k)%Slip(PRN)  ! write O-C residuals
                            write(X38ID,'(I4,F10.1,2I3,2X,A4,2X,I3,A4,F10.3,I2)') Obsweek-1356, ObsSec, 5, IorQ, &
                                    STA%STA(k)%NAME , PRN, "LC", L(N), CycleSlip(k)%Slip(PRN)
                        else
                            
                        end if
                        write(unit=LogID,fmt='(A6,1X,A1,I2,F8.2,3F15.3,E15.7,3F8.3,2F10.3,2X,3F7.3,3I8)') 'PRN',System,PRN_S,Ele, Range, Phase,s,Sat_Clk, &
                                        STD, Ion, rela,  L(N-1), L(N),Res(1)%L(PRN*2),ESC(PRN+56),OrbCorr(PRN+56,1),int(Res(1)%sow),int(PCORSow),int(toe)
                    end if
                end do ! do i=1,ObsData%PRNS

                ! cycle
                if ((STA%STA(k)%SKD=="F") .or. (STA%STA(k)%SKD=="f")) then
                    Res(k)%N=N
                    Res(k)%PRN=PRNPRN
                    Res(k)%Code=Code
                    Res(k)%L=L
                    cycle
                else     ! For the kinematic station
                     Num=N
                     
                    ! ********Eliminate the coordinate  parameter every epoch*******
                    if (Pos_State=="K") then
                        if (ADmethod=='LS') then
                            call Elimi_Para(Nbb,U,SatNum+ParaNum,1)
                            call Elimi_Para(Nbb,U,SatNum+ParaNum,2)
                            call Elimi_Para(Nbb,U,SatNum+ParaNum,3)
                        elseif (ADmethod=='KF') then
                            call KF_Change(InvN,dx,SatNum+ParaNum,1,'pos')
                            call KF_Change(InvN,dx,SatNum+ParaNum,2,'pos')
                            call KF_Change(InvN,dx,SatNum+ParaNum,3,'pos')
                        end if
                        CuParaNum=CuParaNum+3
                    elseif (Pos_State=="S") then
                        if (dabs(InvN(i,i))==0.d0) CuParaNum=CuParaNum+3
                    elseif (Pos_State=="F") then
                        A(1:N,1:3)=0.d0
                    end if
                    
                    ! ********Eliminate the tropsphere parameter *******
                    if ((TropLen/=0.d0) .and. (ADmethod=='LS'))  then
                        if ((Obsweek-TropWeek)*604800.d0+ObsSec-TropSow>=TropLen) then
                            TropWeek=Obsweek
                            TropSow=ObsSec
                            CuParaNum=CuParaNum+1
                            if (Nbb(4,4)==0.d0) then  ! first epoch
                                Nbb(4,4)=1.d4  ! tropsig=0.01m
                                U(4)=1.d4*dx(4)
                            else
                                call Elimi_Para(Nbb,U,SatNum+ParaNum,4)
                                Nbb(4,4)=4.d4
                                U(4)=4.d4*dx(4)
                            end if
                        end if
                    elseif ((TropLen/=0.d0) .and. (ADmethod=='KF'))  then
                        call KF_Change(InvN,dx,SatNum+ParaNum,4,'zwd')
                        if (TropWeek/=0) then  !   http://www.navipedia.net/index_php/Parameters_adjustment_for_PPP
                            InvN(4,4)=InvN(4,4)+1.d-4/3600.d0*((Obsweek-TropWeek)*604800.d0+ObsSec-TropSow )
                        end if
                        TropWeek=Obsweek
                        TropSow=ObsSec
                    end if

                     ! ********Eliminate the rec_clk parameter every epoch*********
                    if (ADmethod=='LS') then
                        call Elimi_Para(Nbb,U,SatNum+ParaNum,5)
                        if ( (ParaNum>5) .and. ((index(isb_mode,"WN") /=0) .or. (index(isb_mode,"wn") /=0)) )then  ! ISB estimated, white noise
                            do i=6, ParaNum
                                call Elimi_Para(Nbb,U,SatNum+ParaNum,i)
                                CuParaNUm=CuParaNUm+1
                            end do
                        end if
                    elseif (ADmethod=='KF') then
                        call KF_Change(InvN,dx,SatNum+ParaNum,5,'clk')
                        if ( (ParaNum>5) .and. ((index(isb_mode,"WN") /=0) .or. (index(isb_mode,"wn") /=0)) )then  ! ISB estimated, white noise
                            do i=6, ParaNum
                                call KF_Change(InvN, dx,SatNum+ParaNum, i, 'clk')
                                CuParaNUm=CuParaNUm+1
                            end do
                        elseif ( (ParaNum>5) .and. ((index(isb_mode,"RW") /=0) .or. (index(isb_mode,"rw") /=0)) ) then     ! Set ISB as a random walk noise, (2m)**2/day
                            do i=6, ParaNum
                                if (GLOIFB==1 .and. i==6+INT_SystemUsed(1)) then ! If linear model of GLONASS IFB
                                    cycle
                                elseif (SystemUsed(1)==.false. .and. GLOIFB==14 .and. i<=19) then ! if frequency depended of GLONASS-only IFB
                                    N=N+1
                                    A(N, 6:ParaNum)=1.d2  ! Tight constraint of Sum(IFB)=0.0
                                    L(N)=0.d0
                                    PP(N)=1.d0
                                    if (InvN(i,i)==0.d0)    InvN(i,i)= 10.d0**2  ! GLONASS IFB frequency depend
                                end if
                                if (dabs(InvN(i,i))>0.d0) then
                                    InvN(i,i)=InvN(i,i) +isb_step**2/86400.d0*Interval
                                end if
                            end do
                        end if
                    end if
                    CuParaNUm=CuParaNUm+1

                    ! ********Random walk of ionosphere parameter*********
                    if (If_Est_Iono .and. (index(ObsCombine,"P1") /=0 .or. index(ObsCombine,"P2") /=0 .or. index(ObsCombine,"P3") /=0)) then
                        do i=1, SatNum
                            if (InvN(ParaNum-SatNum+i,ParaNum-SatNum+i)>0.d0 .and. any(A(1:N, ParaNum-SatNum+i)/=0.d0) ) then
                                InvN(ParaNum-SatNum+i,ParaNum-SatNum+i)=InvN(ParaNum-SatNum+i,ParaNum-SatNum+i)+4.d0/3600.d0*Interval
                            elseif (any(A(1:N, ParaNum-SatNum+i)/=0.d0)) then
                                InvN(ParaNum-SatNum+i,ParaNum-SatNum+i)=0.5d0**2
                            end if
                        end do
                    end if

                     if ((Sta%FixNum>=0) .and. (proc_mod==3)) then
                         call Net_Diff(A(1:N,:), L(1:N),PRNPRN(1:N),Code(1:N),N,ParaNum+SatNum, Num,k, epoch)   ! net difference
                     end if
                end if  !  if ((STA%STA(k)%SKD=="F") .or. (STA                     

                call InvSqrt(MATMUL(Transpose(B(1:NN,1:4)), B(1:NN,1:4)), 4 ,invN2(1:4,1:4))
                PDOP=0.d0
                do i=1,3
                    PDOP=PDOP+invN2(i,i)
                end do
                PDOP=dsqrt(PDOP)
                if (PDOP>MaxPDOP .or. isnan(PDOP)) then
                    write(unit=LogID,fmt='(5X,A5,F5.1,A15)') 'PDOP=', PDOP, '>maxPDOP, skip'
                    cycle
                end if

                do i=1,N
                    write(unit=LogID,fmt='(A10,I3,A3,F10.3)') 'newO-C',PRNPRN(i),Code(i),L(i)
                    L(i)=L(i)*PP(i)
                end do
                
                do i=6,ParaNum
                    if (all(A(1:N,5)-A(1:N,i)==0.d0)) then ! Only one system
                        A(1:N,i)=0.d0
                    end if
                end do
                
                if (ADmethod=='KF') then
                    call InvSqrt(InvN, SatNum+ParaNum ,Nbb)   ! 历史法方程信息
                    U=matmul(Nbb, dx)
                end if
                ! *********** For The Normal Equation ************
                do i=1,N
                    call Change_NEQ(Nbb, U, SatNum+ParaNum, A(i,:), L(i), "add")
                end do


                if (Num<=CuParaNum) then
                    write(unit=LogID,fmt='(5X,A40)') 'Too few observation in this epoch, skip.'
                    cycle
                end if
!                if ( (Num-(CuParaNum-4)*2<4) .and. (epoch>100) ) then
!                    stop   ! 如果至少N-2颗卫星重启，终止程序
!                end if
!                if (epoch<3)  cycle
                Ad_Flag=.true.
                do while (Ad_Flag)
                    Ad_Flag=.false.
!                    Nbb=MATMUL(Transpose(A(1:N,:)), A(1:N,:))
                    call InvSqrt(Nbb, ParaNum+SatNum ,invN)
!                    U=MATMUL(Transpose(A(1:N,:)), L(1:N))
                    dx=MATMUL(invN, U)
                    if (any(isnan(dx))) then
                        write(*,*)  "-------ERROR-------: dx=NAN, please check"
                        write(LogID, '(5X,A40)')  "-------ERROR-------: dx=nan, please check"
                        stop
                    end if
                    V(1:N)=MATMUL(A(1:N,:), dx)-L(1:N)
                    sigma0=dsqrt(DOT_PRODUCT(V(1:N), V(1:N))/(Num-CuParaNum) )
                                   
                    maxV=maxval(dabs(V(1:N)))
                    maxL=maxloc(dabs(V(1:N)),dim=1)
!                    outlier=25.d0/weight
!                    if (outlier<0.10d0) outlier=0.10d0
                    outlier=15.d0
                    if  (Combination(2)==.false.) then
                        if ( (sigma0>2.5d0) .or. ((dabs(maxV)>3.d0*sigma0) .and. (sigma0>1.d0)) ) then    ! 适用于伪距定位
                            goto 125
                        end if
                    elseif ( (dabs(maxV)>outlier)  ) then ! 适用于相位
                        125 call Change_NEQ(Nbb, U, SatNum+ParaNum, A(maxL,:), L(maxL), "cut")
                        A(maxL,:)=0.d0
                        L(maxL)=0.d0
                        Num=Num-1
                        if (Num<CuParaNum) exit
                        Ad_Flag=.true.
                        write(unit=LogID,fmt='(5X,A5,3F10.3,I3,F10.3,A8)') '!!!dx', dx(1:3),PRNPRN(maxL), maxV,'outlier'
                    end if
                    write(unit=LogID,fmt='(5X,A5,3F15.3)') '!!dx',dx(1:3)
                    write(unit=LogID,fmt='(5X,A5,3F15.3,A7,F10.3)') '!!XYZ',AppCoor+dx(1:3),'sigma', sigma0
                end do  ! do while (Ad_Flag)
!                if (Num<CuParaNum) cycle
                do i=1,N
                    if (V(i)/=0.d0) then
                        write(LogID,"(A7,A3,I3,F10.3,F10.2)")  "RES",Code(i),PRNPRN(i), V(i)/PP(i),EleEle(i) ! write observation residuals
                    end if
                end do
                !write(LogID,"(A5, 5F8.3)") "dx:", dx(1:5)
                Coor=STA%STA(k)%Coor(1:3)+dx(1:3)
                if (Pos_State=='K') STA%STA(k)%XYZ=Coor    ! If kinematic, update the approximate coordinate
                call XYZ2BLH(Coor(1), Coor(2), Coor(3), BLH(1), BLH(2), BLH(3))
                Rec_Clk=Rec_Clk+dx(5)/c
                if ((If_ISB))  ISB(1)=dx(6)/c*1.d9
                if ((If_ISB))  ISB(2)=dx(7)/c*1.d9
            
            EpochUsed=EpochUsed+1
            Mean_Coor=Mean_Coor+1.d0/EpochUsed*(Coor-Mean_Coor)  ! Mean coordinate of the station
        
            NEU=MATMUL(Rotation, Coor-STA%STA(k)%TrueCoor)
            UERE=dsqrt(DOT_PRODUCT(NEU,NEU))/PDOP
            write(unit=LogID,fmt='(5X,A5,3F10.3)') '!!NEU',NEU
            write(CoorID,"(I5,I6,4I4,F5.1,4F12.4,I4, 2F8.2,F8.3,I3,3F15.3, 2F15.9,F15.3)") mod(epoch,100000),ObsData%Year, ObsData%Mon, ObsData%Day,  &
     &       ObsData%Hour, ObsData%Min, ObsData%Sec, NEU,dsqrt(DOT_PRODUCT(NEU,NEU)), 0, PDOP, dx(4)+ZHD+ZWD, ObsNum, Coor, BLH
            RMS(1)=RMS(1)+NEU(1)*NEU(1)
            RMS(2)=RMS(2)+NEU(2)*NEU(2)
            RMS(3)=RMS(3)+NEU(3)*NEU(3)
            Mean_NEU=Mean_NEU+1.d0/EpochUsed*(NEU-Mean_NEU)  ! Mean coordinate in NEU direction
        
        
        ! *******End of each station********
    end do  ! do while(flag)
    close(ResO_CID) ! O-C
    close(PcorID) ! O-C
    close(orbcorrID) ! O-C
    end do  ! do k=1,STA%Num
    200 flag=.false.
    
    RMS(1)=dsqrt(RMS(1)/dble(EpochUsed))
    RMS(2)=dsqrt(RMS(2)/dble(EpochUsed))
    RMS(3)=dsqrt(RMS(3)/dble(EpochUsed))

    write(CoorID,"(A)") "--Coordinate====================================================================="
    write(CoorID,"(5X,4F14.4,A20)") RMS,dsqrt(DOT_PRODUCT(RMS,RMS)) ,"//RMS in NEU"
    write(CoorID,"(5X,3F14.4,A35)") Mean_Coor, "//Mean Coordinate"
    write(CoorID,"(5X,3F14.4,A35)") Coor, "//Final Coordinate"
    write(CoorID,"(5X,3F14.4,A35)") dx(1:3), "//Final XYZ Error"
    write(CoorID,"(5X,3F14.4,A35)") NEU, "//Final NEU Error"
    write(CoorID,"(5X,3F14.4,A35)") Mean_Coor - AppCoor,"//Mean XYZ Error"
    write(CoorID,"(5X,4F14.4,A20)") Mean_NEU, dsqrt(DOT_PRODUCT(Mean_NEU,Mean_NEU)),"//Mean NEU Error"
    
    return
    300 write(*,*) "PCOR file file not found, please check."
    pause
    stop
    400 write(*,*) "Orbit correction file not found, please check."
    pause
    stop
    500 write(*,*) "Zone correction file not found, please check."
    pause
    stop
end subroutine