! ================ Process_DD ===================
! PURPOSE
!         Positioning using double difference
! 
! REFERENCE:
!   阳仁贵，袁运斌等. 相位实时差分技术应用于飞行器交会对接研究[J]
!   中国科学：物理学 力学 天文学. 2010,40(5):651~657.
!
! WRITTEN BY: Yize Zhang
! ======================================r=======

subroutine Process_DD
use MOD_FileID
use MOD_Var
use MOD_ObsHead
use MOD_ObsData
use MOD_CycleSlip
use MOD_EOP
use MOD_Rotation
use MOD_STA
use MOD_ZD
use MOD_SD
use MOD_DD
use MOD_NEQ
use MOD_Epo_NEQ
use MOD_NEQ_DP
implicit none
    type(type_ObsData) :: ObsData(2) 
    type(type_ZD) :: ZD(2)
    type(type_SD) :: SD
    type(type_DD) :: DD
    type(type_NEQ) :: NEQ
    type(type_Epo_NEQ) :: Epo_NEQ
    ! CRS to TRS
    real(8) :: MJD,kmjd, Xp,Yp,dUT1,DX00,DY00
    real(8) :: TT
    logical :: flag
    integer :: epoch,Obsweek
    real(8) :: Obssec
    integer :: error, data_flag
    integer(1) :: Flag_Sln
    integer :: i, j, RefSat(5), sys
    integer :: EpochUsed
    real(8) :: dx_Coor(4), Coor(3), BLH(3), Mean_Coor(3), NEU(3), RMS(3), Mean_NEU(3)
    real(8) :: PDOP, Nbb(3,3), InvN(3,3)
    real(8) :: Vel(3)

    ParaNum=3  ! no tropsphere and clock for single system double difference
    if (If_Est_Iono) then       ! If estimate ionosphere
        IonoNum=MaxPRN
        AmbID(1)=FileID_Mark
        FileID_Mark=FileID_Mark+1
    end if
!    open(unit=AmbID(1), file='E:\TUMSAT\long_baseline3\dztd_koma_143',action="read")
!    AmbID(2)=FileID_Mark
!    FileID_Mark=FileID_Mark+1
!    open(unit=AmbID(2), file='E:\TUMSAT\long_baseline3\dztd_kaiy_143',action="read")
    if (TropLen/=0.d0) then       ! If estimate troposphere
        ParaNum=ParaNum+1
    end if
    allocate(ZD(1)%A(MaxPRN,ParaNum), ZD(2)%A(MaxPRN,ParaNum), SD%A(MaxPRN,ParaNum), DD%A(MaxPRN,ParaNum) )
    allocate(NEQ%Ap1(MaxPRN,ParaNum),NEQ%Ap2(MaxPRN,ParaNum),NEQ%Awl(MaxPRN,MaxPRN*2+ParaNum),NEQ%Aw4(MaxPRN,MaxPRN*2+ParaNum))
    allocate(NEQ%Aewl(MaxPRN,ParaNum))
    allocate(NEQ%Nbb(MaxPRN*2+ParaNum,MaxPRN*2+ParaNum),NEQ%InvN(MaxPRN*2+ParaNum,MaxPRN*2+ParaNum))
    allocate(NEQ%U(MaxPRN*2+ParaNum),NEQ%dx(MaxPRN*2+ParaNum) )
    NEQ%Ap1=0.d0; NEQ%Lp1=0.d0; NEQ%Ap2=0.d0; NEQ%Lp2=0.d0
    NEQ%Awl=0.d0; NEQ%Lwl=0.d0; NEQ%Aw4=0.d0; NEQ%Lw4=0.d0
    NEQ%Nbb=0.d0; NEQ%InvN=0.d0; NEQ%U=0.d0; NEQ%dx=0.d0; NEQ%N=MaxPRN*2+ParaNum
    allocate(Epo_NEQ%Ap1(MaxPRN,ParaNum+3*IonoNum),Epo_NEQ%Ap2(MaxPRN,ParaNum+3*IonoNum),Epo_NEQ%Al1(MaxPRN,ParaNum+3*IonoNum),Epo_NEQ%Al2(MaxPRN,ParaNum+3*IonoNum))
    allocate(Epo_NEQ%Awl(MaxPRN,ParaNum+3*IonoNum),Epo_NEQ%Aw4(MaxPRN,ParaNum+3*IonoNum))
    allocate(Epo_NEQ%Nbb(ParaNum+3*IonoNum,ParaNum+3*IonoNum),Epo_NEQ%InvN(ParaNum+3*IonoNum,ParaNum+3*IonoNum),Epo_NEQ%U(ParaNum+3*IonoNum),Epo_NEQ%dx(ParaNum+3*IonoNum) )
    Epo_NEQ%Ap1=0.d0; Epo_NEQ%Lp1=0.d0; Epo_NEQ%Ap2=0.d0; Epo_NEQ%Lp2=0.d0
    Epo_NEQ%Al1=0.d0; Epo_NEQ%Ll1=0.d0; Epo_NEQ%Al2=0.d0; Epo_NEQ%Ll2=0.d0
    Epo_NEQ%Awl=0.d0; Epo_NEQ%Lwl=0.d0; Epo_NEQ%Aw4=0.d0; Epo_NEQ%Lw4=0.d0
    Epo_NEQ%Nbb=0.d0; Epo_NEQ%InvN=0.d0; Epo_NEQ%U=0.d0; Epo_NEQ%dx=0.d0; Epo_NEQ%N=ParaNum+3*IonoNum


    flag=.true.
    epoch=0
    RefSat=0
    EpochUsed=0
    RMS=0.d0
    Mean_Coor=0.d0
    Mean_NEU=0.d0

    ! ***************** Loop epoch-by-epoch ******************
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
        else
            Obssec=Obssec+Interval
            if (Obssec>=604800.0d0) then
                Obsweek=Obsweek+1
                Obssec=Obssec - 604800.d0
            end if
        end if
        if ( (Obsweek-GPSweek_st)*604800.0d0+Obssec-GPSsec_st<0.d0 ) cycle
        if ( (Obsweek-GPSweek_end)*604800.0d0+Obssec-GPSsec_end>0.d0 ) exit
        
        if (If_Est_Iono) then  ! For long baseline RTK, should consider solid tide and ocean load correction
            ! Calculate the transfor matrix between TRS and CRS
            MJD=Obsweek*7.d0+44244.0d0+Obssec/86400.d0 !-Leap_sec/86400.d0 ! In GPST
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
        end if
        
        ! ============= Read obs data at the given time ================
        Data_flag=0
        do i=1,STA%Num  ! STA%Num=2 in double difference
            call ReadObsData(Obsweek,Obssec,ObsData(i),i)
            read(unit=ObsID(i),fmt="(2X)",iostat=error)
            backspace(ObsID(i))
            if (error /=0) flag=.false.  ! Check that if reach the end of observation file.
            if (i==1 .and. dabs((Obsweek-ObsData(i)%Week)*604800.d0+Obssec-ObsData(i)%Sow)<Interval) then
                ObsData(i).Flag=0 ! If reference staiton and time difference less than Interval
            end if
            if (ObsData(i).Flag /=0) then
                Data_flag=1
                ObsData(i).PRNS=0
            end if
        end do
        if (Data_flag==1) cycle
        
        if (mod(epoch,10)==1) write(*,'(A5,I6,2X,I4,4I3,F5.1)') "Epoch", epoch, ObsData(2)%year,ObsData(2)%mon,&
                       ObsData(2)%day, ObsData(2)%hour, ObsData(2)%min, ObsData(2)%sec
        if (ObsHead(2)%Version==2) then
            write(unit=LogID,fmt='(A5,I5,2X, I4,4I3, F5.1)') 'Epoch', epoch, ObsData(2)%year,ObsData(2)%mon,&
                           ObsData(2)%day, ObsData(2)%hour, ObsData(2)%min, ObsData(2)%sec
        else
            write(unit=LogID,fmt='(A5,I5,2X, I4,4I3.2, F5.1)') 'Epoch', epoch, ObsData(2)%year,ObsData(2)%mon,&
                           ObsData(2)%day, ObsData(2)%hour, ObsData(2)%min, ObsData(2)%sec
        end if
        write(unit=CSID,fmt='(A5,I5,2X, I4,4I3.2, F5.1)') 'Epoch', epoch, ObsData(2)%year,ObsData(2)%mon,&
                       ObsData(2)%day, ObsData(2)%hour, ObsData(2)%min, ObsData(2)%sec

        ! ===== For each station, form zero difference ========
        if (Pos_State=='K') STA%STA(2)%Coor=Coor  ! If kinematic, update the approximate coordinate
        do i=1,STA%Num  ! STA%Num=2 in double difference
!            if (any(STA%STA(i)%TrueCoor==0.d0) .or. (Pos_State=='K')) then
!            if (any(STA%STA(i)%TrueCoor==0.d0) .or. ((Pos_State=='K') .and. .not.(STA%STA(i)%flag_InitialCoor)) ) then
            if ( any(STA%STA(i)%TrueCoor==0.d0) .or. (Pos_State=='K' .and. any(dabs(STA%STA(i)%TrueCoor-STA%STA(i)%Coor)>50.d0))  ) then
!            if (any(STA%STA(i)%TrueCoor==0.d0)) then
                if (.not. Combination(3)) then  ! If no doppler velocity
                    call Bancroft(ObsData(i), i, STA%STA(i)%Coor)   ! If Kinematic
                    if (any(STA%STA(i)%Coor==0.d0)) cycle
                elseif (i>1 .and. Combination(3)) then
!                    call polyfit(NEQ_DP%Sow(2:5)-ObsSec, NEQ_DP%Vel(:,2:5), 4, NEQ_DP%Flag_Sln(2:5), 0.d0, Vel, j)  ! Predict velocity by Quadratic polyfit doppler velocity
!                    if (all(Vel/=0.d0)) then
!                        NEQ_DP%dt=Obssec-NEQ_DP%Sow(j)
!                        STA%STA(i)%Coor=NEQ_DP%Coor+Vel*NEQ_DP%dt
                    NEQ_DP%dt=Obssec-NEQ_DP%Sow(5)
                    if (all(NEQ_DP%Vel(:,5)/=0.d0) .and. NEQ_DP%dt<=5.d0) then  ! If Doppler velocity is valid and ambiguity fixed, use it for the new float position
                        STA%STA(i)%Coor=NEQ_DP%Coor+NEQ_DP%Vel(:,5)*NEQ_DP%dt
                        Vel_Used=1
                    else
                        call Bancroft(ObsData(i), i, STA%STA(i)%Coor)   ! If Kinematic
                        if (any(STA%STA(i)%Coor==0.d0)) cycle
                        Vel_Used=0
                    end if
                end if
            end if
            call Zero_Difference(ObsData(i), i, ZD(i))
        end do

        ! ========= Form station difference ============
        call Station_Difference(ZD, SD)

        ! ========== Check reference satellite in this epoch ============
        if ( (arc_epochs/=0) .and. (mod(epoch,arc_epochs)==1) )  DD%RefSat=0  ! 
        if (ar_mode==2)  DD%RefSat=0  ! Instantaneous AR
        call Get_RefSat(SD, DD%RefSat, RefSat)
        if (all(RefSat==0)) then
            write(LogID,'(5X,A30)') '未找到参考星，不计算此历元。'
            cycle          ! If no reference satellite found(almost not possible)
        end if
        
        ! ========= Form double difference =============
        call Double_Difference(SD, DD, RefSat)
        
        if ( (any(DD%RefSat-RefSat/=0)) .and. (DD%PRNS<3+1) ) then ! 如果卫星数太少
            write(LogID,'(5X,A40)') 'Too few sat & no previous ref sat. Skip!'
            cycle
        end if

        ! If new reference satellite found, reset the NEQ
        if (ar_mode/=2) then
            call Reset_NEQ(RefSat, DD, NEQ, Epo_NEQ)
        end if

        if (.not.(If_Est_Iono) .and. If_IonoCompensate) then
            call IonoCompensate(Obsweek,Obssec, RefSat, DD,Epo_NEQ) ! Compensate DD ionosphere residuals
        end if
        do sys=1, 5   ! mark the reference satellite
            if (RefSat(sys)/=0) DD%RefSat(sys)=RefSat(sys)
        end do
        
        ! ======== Form normal equation ===============
        call Form_NEQ(SD, DD, NEQ, Epo_NEQ)

        Nbb = matmul(  matmul( transpose(NEQ%Ap1(1:NEQ%PRNS,1:3)), NEQ%P(1:NEQ%PRNS,1:NEQ%PRNS) ), NEQ%Ap1(1:NEQ%PRNS,1:3)  )
        call Invsqrt(Nbb, 3, InvN)
        PDOP=0.d0
        do i=1,3
            PDOP=PDOP+InvN(i,i)
        end do
        PDOP=dsqrt(PDOP)/100.d0

        !  ======== Least squart estimation for the coordinate ========
        if (DD%PRNS<4) then
            write(unit=LogID,fmt='(A10,A30)') '%%NEU','Insufficient Sat Num '
            cycle
        end if
        if (If_Est_Iono .and. IonoNum>0) then
            call Solve_NEQ_Iono(NEQ, Epo_NEQ, dx_Coor, Flag_Sln)
        else
            call Solve_NEQ(NEQ, Epo_NEQ, dx_Coor(1:3), Flag_Sln)
        end if
        if (.not.(If_Est_Iono) .and. If_IonoCompensate .and. Flag_Sln<3) then
            call IonoUpdate(ObsWeek, ObsSec, NEQ, Epo_NEQ) ! Update DD ionosphere residuals
        end if

        ! ================= Write the result ===============
        Coor=STA%STA(2)%Coor+dx_Coor(1:3)
        call XYZ2BLH(Coor(1), Coor(2), Coor(3), BLH(1), BLH(2), BLH(3))
        if (Flag_Sln<=3) then
            NEQ_DP%Coor=Coor
            NEQ_DP%Vel=eoshift(NEQ_DP%Vel,shift=1,dim=2, boundary=NEQ_DP%dx(1:3))
            NEQ_DP%Sow=eoshift(NEQ_DP%Sow,shift=1,boundary=ObsSec)
            NEQ_DP%Flag_Sln=eoshift(NEQ_DP%Flag_Sln,shift=1, boundary=Flag_Sln)
        end if
        NEU=matmul(STA%STA(2)%Rotation, Coor-STA%STA(2)%TrueCoor)
       
        write(CoorID,"(I5,I6,4I4,F5.1,4F12.4,I4,F8.2,F8.3,I3,3F15.3, 2F15.9,F15.6)") mod(epoch,100000),ObsData(2)%Year, ObsData(2)%Mon, & 
                 ObsData(2)%Day,  ObsData(2)%Hour, ObsData(2)%Min, ObsData(2)%Sec, NEU, dsqrt(DOT_PRODUCT(NEU,NEU)),&
                  Flag_Sln, PDOP, STA%STA(2)%Trop%ZHD+STA%STA(2)%Trop%ZWD+dx_Coor(4),  Epo_NEQ%PRNS, Coor, BLH
         write(unit=LogID,fmt='(A10,3F10.3)') '%%NEU',NEU

        RMS(1)=RMS(1)+NEU(1)*NEU(1)
        RMS(2)=RMS(2)+NEU(2)*NEU(2)
        RMS(3)=RMS(3)+NEU(3)*NEU(3)
        EpochUsed=EpochUsed+1
        Mean_Coor=Mean_Coor+1.d0/EpochUsed*(Coor-Mean_Coor)  ! Mean coordinate of the station
        Mean_NEU=Mean_NEU+1.d0/EpochUsed*(NEU-Mean_NEU)
    end do  ! do while(flag)
    ! ************ End of epoch loop ****************

    RMS(1)=dsqrt(RMS(1)/dble(EpochUsed))
    RMS(2)=dsqrt(RMS(2)/dble(EpochUsed))
    RMS(3)=dsqrt(RMS(3)/dble(EpochUsed))
    write(CoorID,"(A)") "--Coordinate==================================================================="
    write(CoorID,"(5X,4F14.4,A20)") RMS,dsqrt(DOT_PRODUCT(RMS,RMS)) ,"//RMS in NEU"
    write(CoorID,"(5X,4F14.4,A20)") Mean_NEU, dsqrt(DOT_PRODUCT(Mean_NEU,Mean_NEU)),"//Mean NEU Error"
    write(CoorID,"(5X,3F14.4,A35)") Mean_Coor , "//Mean Coordinate"
    write(CoorID,"(5X,3F14.4,A20)") STA%STA(2)%Coor(1:3)+dx_Coor(1:3), "//Final Coor"
    !write(CoorID,"(A25,I6,I6, F8.4)") 'amb fix success num:  ', amb_success, Epochused, dble(amb_success/Epochused)

end subroutine
