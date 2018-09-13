! *********************************************
! Entrance of the program Net_Diff
!
! NOTICE: 
!           This program is for static or kinematic SPP/PPP/RTK for 
!     single-station. As for multi-station PPP, net-combination
!     PPP or satellite clock estimation, please turn to Net_PPP.
!            In the other hand, this programe can used for SPP/PPP
!    with augmentation information or station/double difference
!   solution.
!
!  Any discussion or question, please contact with me.
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ******************** ^ ** ^ ******************
! ********************** V ********************


program main
use MOD_FileID
use MOD_FileDir
use MOD_STA
use MOD_Time
Use MOD_ObsHead
use MOD_Var
use MOD_EOP
use MOD_Constant
use MOD_Ant
implicit none
    character(200) :: ObsFile, NavFile
    type(cal_time) t
    character(7) :: str_day
    character(4) :: sta_name
    integer ::  i
    integer :: year, doy, GPSweek, GPSday, int_day
    integer :: mon,day,hour,min
    real(8) :: sec
    character(5) :: str_GPST
    character(2) :: temp
    character(1) :: temp2
     real(8) :: MJD
     logical :: alive

    call Control_Net
    call Coor_Table

    int_day=int_year*1000+int_doy
    write(str_day,"(I7)") int_day

    call Day2GPST(int_year, int_doy, GPSweek, GPSday)    ! GPS week and GPS day
    call UTC2MJD(int_year,1,int_doy, 12, 0, 0.d0, MJD)   ! In the middle of the day
    write(str_GPST,"(I4,I1)") GPSweek,GPSday    
    
    call CPU_time(t%t_begin) 
    
    ! File preparation
    if  ( (index(Orbit,"BRD")/=0) .or. (index(Clk,"BRD")/=0) .or. ((index(ObsCombine,"PC")==0) .and. (iontype==1)) ) then  ! If broadcast navigation
        if (SystemUsed(1) .or. SystemUsed(5) .or. SystemUsed(6) ) then  ! GPS or QZSS or IRNSS
            NavFile=trim(NavDir)//"brdm"//str_day(5:7)//"0."//str_day(3:4)//"p"
            call ReadNavHead(NavFile)
            call ReadNavData
            if (IsCNAV) then  ! Read CNAV ephemeris of GPS
                NavFile=trim(NavDir)//"brdx"//str_day(5:7)//"0."//str_day(3:4)//"x"
                call ReadNavHead(NavFile)
                call ReadNavData
            end if
        end if
        if (SystemUsed(2)) then
            NavFile=trim(NavDir)//"brdm"//str_day(5:7)//"0."//str_day(3:4)//"p"
            call ReadNavHead_R(NavFile)
            call ReadNavData_R
        end if
        if (SystemUsed(3)) then
            if (IorQ==0) then
                NavFile=trim(NavDir)//"brdm"//str_day(5:7)//"0."//str_day(3:4)//"p"
            elseif (IorQ==3) then
                NavFile=trim(NavDir)//"brdm"//str_day(5:7)//"0."//str_day(3:4)//"pQ"
            end if
            call ReadNavHead_C(NavFile)
            call ReadNavData_C
            call Get_BEB
!            NavFile="D:\Program\RDRN2RNX\SatClkPara"
!            call ReadNavClk(NavFile)
        end if
        if (SystemUsed(4)) then
            NavFile=trim(NavDir)//"brdm"//str_day(5:7)//"0."//str_day(3:4)//"p"
            call ReadNavHead_E(NavFile)
            call ReadNavData_E
        end if
    end if
    if  (index(Orbit,"SP3")/=0) then  ! If precise orbit
        SP3File=trim(sp3Dir)//str_GPST//".sp3"
         ! Read SP3Head and 10 epochs of SP3 data
        call ReadSP3Head(SP3File)
        call ReadSP3Data(10)
!        do i=1,14
!            write(temp,'(I2.2)') i
!            SP3File=trim(sp3Dir)//"11600_S"//temp//".rpt15"
!            inquire(file=SP3File,exist=alive)
!            if (.not. alive) then
!                write(*,*) "SP3 flie: """//trim(SP3File)//""" doesn't exist!"
!            end if
!            open(unit=100+i, file=SP3File,action='read')
!        end do
!        SP3Time= -14.d0
!        do i=0,9
!            call ReadMultiOrb(1, 1841, 348300.d0+i*300.d0)
!        end do
    end if
    if (index(Clk,"CLK")/=0)  then  ! If precise clock
        ClkFile=trim(clkDir)//str_GPST//".clk"
        ! Read clock file header and two epochs of clk data
        call ReadClkHead(ClkFile)
        call ReadClkData
        call ReadClkData
    end if
     ! Read ObsHead
     if ( (ObsType=='X71') .or. (ObsType=='x71') ) then
         do i=1,STA%Num
             ObsFile= trim(ObsDir)//'x71_'//STA%STA(i)%Name(1:2)//"_"//STA%STA(i)%Name(3:4)//"@"//str_day
             ! inquire that if the obs file exists
             inquire(file=ObsFile,exist=alive)
             if (.not. alive) then
                 write(*,*) "Observation flie: """//trim(ObsFile)//""" doesn't exist!"
                 pause
                 stop
             end if
             ObsID(i)=FileID_Mark
             FileID_Mark=FileID_Mark+1
             open(unit=ObsID(i), file=ObsFile,action="read")
         end do
         TimeSys="BDT"
     elseif ( (ObsType=='X11') .or. (ObsType=='x11') ) then
        do i=1,STA%Num
             ObsFile= trim(ObsDir)//'gps_'//STA%STA(i)%Name(1:2)//"_"//STA%STA(i)%Name(3:4)//"@"//str_day
             ! inquire that if the obs file exists
             inquire(file=ObsFile,exist=alive)
             if (.not. alive) then
                 write(*,*) "Observation flie: """//trim(ObsFile)//""" doesn't exist!"
                 pause
                 stop
             end if
             ObsID(i)=FileID_Mark
             FileID_Mark=FileID_Mark+1
             open(unit=ObsID(i), file=ObsFile,action="read")
         end do
        TimeSys="BDT"
     else
        do i=1,STA%Num
            ObsFile= trim(ObsDir)//STA%STA(i)%Name//str_day(5:7)//"0."//str_day(3:4)//"o"
            call ReadObsHead(ObsFile,i)
            STA%STA(i)%Ant_Type=ObsHead(i)%Ant_Type
        end do
    end if
    if (TimeSys=="GPS") then
        ObsTime=0.d0
    elseif ((TimeSys=="BDS") .or. (TimeSys=="BDT")) then
        ObsTime= 14.d0
    else
        write(*,*) "Obsveration time system unknown, will default as GPS."
!        pause
        ObsTime=0.d0
    end if

    ! Get DCB
    call Get_DCB
    
    ! Ionospheric Delay
    call GPST2UTC(GPSweek, GPSday*86400.d0+0.1d0, MJD,year,mon,day,hour,min,sec)
    write(str_ymd,'(I4,I2.2,I2.2)') year, mon, day
    if ( (index(ObsCombine,"PC")==0).and. (index(ObsCombine,"LC")==0) )then
        if (iontype==2) then ! GIM model
            IONFile=trim(IonDir)//str_day(5:7)//'0.'//str_day(3:4)//'i'
            call ReadIONEX(IONFile)
        elseif (iontype==3) then ! Get 14 Parameter of BeiDou Ionospheric Delay
            IONFile=trim(IonDir)//str_ymd//'_14Psim.dat'
            call Get_Iono14(IONFile)
        elseif (iontype==4) then ! grid model for BeiDou
             call rd_igpdely_online(year,mon,day)
        end if
    end if

     ! Get EOP parameter
    call UTC2MJD(int_year,1,int_doy,12, 0, 0.d0, MJD)   ! In the middle of the day
    call Get_EOP(MJD, EOP%MJD, EOP%X,EOP%Y,EOP%dUT1, EOP%dX, EOP%dY)
    
     ! Get antenna corrections information
    call Get_Ant(MJD)
    if (JNum>0 .and. index(sp3dir(len_trim(sp3dir)-2:len_trim(sp3dir)),'gbm') /=0) then
!        Ant(2+GNum+RNum+CNum+NumE)%PCO(1:3,1)=4.16d0+3.2d0 !??? absolute value is wrong
!        Ant(2+GNum+RNum+CNum+NumE)%PCO(1:3,2)=4.16d0+3.2d0 !???
!        Ant(3+GNum+RNum+CNum+NumE)%PCO(1:3,1)=4.7d0+3.2d0
!        Ant(3+GNum+RNum+CNum+NumE)%PCO(1:3,2)=4.7d0+3.2d0
    end if
    if (CNum>0 .and. index(sp3dir(len_trim(sp3dir)-2:len_trim(sp3dir)),'gbm') /=0 .and. int_year*1000+int_doy>2014197) then
        call Get_BDS_PCO_gbmwum(int_year,int_doy, 1)
    elseif (CNum>0 .and. index(sp3dir(len_trim(sp3dir)-2:len_trim(sp3dir)),'wum') /=0 ) then
        call Get_BDS_PCO_gbmwum(int_year,int_doy, 2)
    end if

    ! Get GLONASS frequency and leap second
    call Get_GLO_Fre 
    
    ! Get Ocean load coefficient. If coordinate is 0.000, it will be 0.000. But this is usually in kinematic mode, it is OK we ignore ocean load correction.
    ! Or we can give a initial coordinate
    do i=1,STA%Num
        call Get_OceanLoad_Coe(STA%STA(i)%BLH(1), STA%STA(i)%BLH(2), STA%STA(i)%Name, STA%STA(i)%OLC%CMC, STA%STA(i)%OLC%OLC,STA%STA(i)%OLC%found)
    end do

    CoorID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    write(temp,"(I1)") Sta%FixNum
    write(temp2,'(I1)') int(delay/60)
    open(unit=CoorID, file=trim(OutDir)//"Coor_"//str_day//"_"//STA%STA(1)%Name//"-"//STA%STA(STA%Num)%Name//".txt",action="write",err=100)
    write(CoorID,"(A40)") "===============NET Diff================="
    write(CoorID,"(A40)") "Developed by Yize Zhang, zhyize@163.com"
    write(CoorID,"(A15,A7)") "day: ",str_day
    write(CoorID,"(A20)") "Station Fixed: "
    do i=1,STA%FixNum
        write(CoorID,"(5X,A5,3F15.3,3F10.3,A25)") STA%STA(i)%Name,STA%STA(i)%TrueCoor, STA%STA(i)%NEU, STA%STA(i)%Ant_Type
    end do
    write(CoorID,"(A20)") "Station Estimated:"
    do i=STA%FixNum+1,STA%Num
        write(CoorID,"(5X,A5,3F15.3,3F10.3,A25)") STA%STA(i)%Name,STA%STA(i)%Coor,STA%STA(i)%NEU, STA%STA(i)%Ant_Type
    end do
    write(CoorID,"(A12,A)") 'eopfile:  ',trim(EOPFile)
    write(CoorID,"(A12,A)") 'antfile:  ',trim(AntFile)
    write(CoorID,"(A12,A)") 'planetfile:  ',trim(PlanetFile)
    !write(CoorID,"(A12,A)") 'tidefile:  ',trim(TideFile)
    write(CoorID,"(A12,A)") 'glofrefile:  ',trim(GloFreFile)
    write(CoorID,"(A12,A)") 'p1c1file:  ',trim(P1C1File)
    write(CoorID,"(A12,A)") 'dcbfile:  ',trim(DCBFile)
    write(CoorID,"(A12,A)") 'obsdir:  ',trim(ObsDir)
    write(CoorID,"(A12,5X,A5)") 'obstype:  ',ObsType
    if (ObsType=='x71') then
        write(CoorID,"(A12,I5)") 'IORQ:  ',IorQ
    end if
    if  (index(Orbit,"BRD")/=0) then
        write(CoorID,"(A12,A)") 'navdir:  ',trim(NavDir)
    elseif (index(Orbit,"SP3")/=0) then
        write(CoorID,"(A12,A)") 'sp3file:  ',trim(SP3File)
    end if
    if (index(Clk,"CLK")/=0) then
        write(CoorID,"(A12,A)") 'clkfile:  ',trim(ClkFile)
    end if
    write(CoorID,"(A12,5X,F5.1)") 'interva(s):  ', interval
    write(CoorID,"(A12,5X,I5)") 'limele:  ', int(limele)
    write(CoorID,"(A12,5X,I5)") 'limSNR:  ', int(limSNR)
    write(CoorID,"(A12,5X,I5)") 'maxpdop:  ', int(MaxPDOP)
    write(CoorID,"(A12,5X,A8,A5,A5)") 'trop_model:  ', trim(cdattype), trim(cztd), trim(cmap)
    if (TropLen>0) then
        write(CoorID,"(A12,5X,I5)") 'trop_len:  ', int(TropLen)
    end if
    call GPST2UTC(GPSweek_st, GPSsec_st, MJD,year,mon,day,hour,min,sec)
    write(CoorID,"(A12,5X,I5,4I3,F5.1,I5,I8)") 't_start: ', year,mon,day,hour,min,sec, GPSweek_st, int(GPSsec_st)
    call GPST2UTC(GPSweek_end, GPSsec_end, MJD,year,mon,day,hour,min,sec)
    write(CoorID,"(A12,5X,I5,4I3,F5.1,I5,I8)") 't_end:  ', year,mon,day,hour,min,sec, GPSweek_end, int(GPSsec_end)
    write(CoorID,"(A12,5X,A5)") 'freq_comb:  ',freq_comb
    if (Combination(1)) then
        write(CoorID,"(A12,5X,A5,F5.1)",advance='no') 'combination:  ','PC',sigPC
    end if
    if (Combination(2)) then
        write(CoorID,"(A5,F7.3)",advance='no') 'LC',sigLC
    end if
    if (Combination(3)) then
        write(CoorID,"(A5,F5.2)",advance='no') 'DP',sigDP
    end if
    write(CoorID,"(A)") ''
    write(CoorID,"(A12,5X,6L5)") 'systemused:  ', SystemUsed(1:6)
    write(CoorID,"(A12,5X,A5)") 'admethod:  ', ADmethod
    write(CoorID,"(A12,5X,A5)") 'pos_state:  ', Pos_State
    !write(CoorID,"(A12,5X,A5)") 'posmethod:  ', PosMethod
    write(CoorID,"(A12,5X,I5)") 'proc_mod:  ', proc_mod
    if ((proc_mod>=1) .and. (proc_mod<=3)) then
        write(CoorID,"(A12,5X,A)") 'pcorfile:  ', trim(pcorfile)
    end if
    if ((proc_mod==2) .or. (proc_mod==3)) then
        write(CoorID,"(A12,5X,A)") 'orbcorrfile:  ', trim(orbcorrfile)
    end if
    if (proc_mod==3) then
        write(CoorID,"(A12,5X,A)") 'zonecorrfile:  ', trim(zonecorrfile)
    end if
    if (proc_mod==5) then
        write(CoorID,"(A12,5X,4F5.1)") 'dd_coe:  ', a1, a2, b1, b2
        if (ar_mode==0) then
            write(CoorID,"(A12,5X,A15)") 'ar_mode:  ', 'Float'
        elseif (ar_mode==1) then
            write(CoorID,"(A12,5X,A15)") 'ar_mode:  ', 'Continuous'
            write(CoorID,"(A12,5X,2I4)") 'fixele:  ', int(FixEle)
        elseif (ar_mode==2) then
            write(CoorID,"(A12,5X,A15)") 'ar_mode:  ', 'Instantaneous'
            write(CoorID,"(A12,5X,2I4)") 'fixele:  ', int(FixEle)
        elseif (ar_mode==3) then
            write(CoorID,"(A12,5X,A15)") 'ar_mode:  ', 'Fix and hold'
            write(CoorID,"(A12,5X,2I4)") 'fixholdele:  ', int(FixEle), int(HoldEle)
        end if
        write(CoorID,"(A12,5X,L5,I3)") 'partial_ar:  ', partial_AR, parARnum
        write(CoorID,"(A12,5X,F5.1)") 'minratio:  ', minratio
        write(CoorID,"(A12,5X,L5)") 'est_wl_amb:  ', If_Est_WL
    end if
    if ((proc_mod>=1) .and. (proc_mod<=3)) then
        write(CoorID,"(A12,5X,I5)") 'delay:  ', int(delay0)
        write(CoorID,"(A12,5X,I5)") 'clktype:  ', clktype
    end if
    write(CoorID,"(A12,5X,L5)") 'isb:  ', If_ISB
    write(CoorID,"(A12,5X,L5)") 'ifb:  ', If_IFB
    write(CoorID,"(A12,5X,A)") 'cycleslip:  ', CSmethod
    write(CoorID,"(A12,5X,A5)") 'obscombine:  ', trim(ObsCombine)
    write(CoorID,"(A12,5X,L5)") 'est_iono:  ', If_Est_Iono
    if (index(ObsCombine,'PC')==0) then
        if (iontype>1) then 
            write(CoorID,"(A12,5X,A)") 'ionfile:  ', trim(IONFile)
        end if
        write(CoorID,"(A12,5X,I5)") 'iontype:  ', iontype
    end if
    write(CoorID,"(A12,5X,A5,A6,A5,I5)") 'smooth:  ', Var_smooth,Smooth_Method, Smooth_Combine,int(Smooth_Time)
    write(CoorID,"(A12,5X,A)") 'outdir:  ', trim(OutDir)
    write(CoorID,"(A)") "--End of header========================================"
    write(CoorID,"(A5,A27,4A12, A4,2A8,A3,6A15)") "Epoch","GPST(yyyymmddhhmmss)","N Error/m","E Error/m","U Error/m","Pos Error/m", &
                           "Q", "PDOP","ZTD","N","X/m(ECEF)","Y/m(ECEF)","Z/m(ECEF)","Lat/deg","Lon/deg","Hgt/m"
     write(CoorID,"(A)") "++Coordinate================================================================="

    LogID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=LogID, file=trim(OutDir)//"Log_"//str_day//"_"//STA%STA(1)%Name//"-"//STA%STA(STA%Num)%Name//".txt",action="write",err=100)
    write(LogID,"(A40)")    "============ Begin Computing ==========="
    
    CSID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=CSID, file=trim(OutDir)//"cs_"// str_day // ".txt",action="write",err=100)

    if ( (proc_mod==1) .or. (proc_mod==2) .or. (proc_mod==3)) then
        ResO_CID=FileID_Mark   ! Correction error, used for smooth comparation
        FileID_Mark=FileID_Mark+1
        open(unit=ResO_CID, file=trim(OutDir)//"Resi_O_C.txt",action="readwrite",err=100)
    
        X38ID=FileID_Mark   ! Correction error, used for smooth comparation
        FileID_Mark=FileID_Mark+1
        open(unit=X38ID, file=trim(OutDir)//"X38_"// str_day // ".txt",action="write",err=100)
    elseif (proc_mod==5) then    
        LambdaID=FileID_Mark
        FileID_Mark=FileID_Mark+1
        open(unit=LambdaID, file=trim(OutDir)//"lambda_"// str_day // ".txt",action="write",err=100)
        write(LambdaID,"(A19)")    "Lambda processing."
    end if

    ! =====Main process=====
    if ( (Var_smooth=="y") .or. (Var_smooth=="Y") ) then
    if  (index(Smooth_Method,"SMT") /=0)  then
        call PreProcess
    end if
    end if

!    call Process_BRD_Sp3  ! navigation
!      call Process_Net
    if (proc_mod==0 .or. proc_mod==6) then  
        call Process  ! SPP(0)/PPP(0)/DSPP(6)/DPPP(6)
    elseif ( (proc_mod==1) .or. (proc_mod==2) .or. (proc_mod==3)) then
        call Process_Corr    ! 加分区改正数或等效钟差
!    elseif (proc_mod==4) then
!        call Process2      ! 2min 站间差
    elseif (proc_mod==5) then
        call Process_DD    ! RTK
    else
        write(*,*) "proc_mod error"
        pause 
        stop
    end if
    ! =======================

    do i=1,STA%Num
        close(ObsID(i))
    end do
    close(LogID)
    close(CSID,status='delete')
    if ( (proc_mod==1) .or. (proc_mod==2) .or. (proc_mod==3)) then
        close(ResO_CID,status='delete') ! O-C
        close(X38ID,status='delete') ! X38
    elseif (proc_mod==5) then
        close(LambdaID, status='delete')
    end if
    close(CoorID)
    
    call CPU_time(t%t_end)!(time_end1)
    write(*,"(A16F10.2,A1)") "Total time:" , t%t_end - t%t_begin,"s"
    stop
    100 write(*,*) '-------ERROR------: Opening output file error, please check if the file exist or it can be read (If Net_Diff.dll exist, kill it)'
    stop

end program