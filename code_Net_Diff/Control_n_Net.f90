! ===========Control_Net============
! PURPOSE:
!      Read the control file and get the 
!      input file directories and other variables!
! 
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ============End of Header============

subroutine Control_Net
use MOD_FileDir
use MOD_FileID
use MOD_Var
use MOD_NavData
use MOD_NavData_R
implicit none
    logical :: alive
    character(300) :: line
    character(100) :: temp
    character(100) :: str_pos_state, str_SystemUsed, str_Combination 
    character(3) :: str_partial_ar, str_ISB,str_IFB, str_IonoCompensate
    integer :: year, mon, day, hour, min, n,PRN, idx
    real(8) :: sec
    
    
    inquire(file=ConFile,exist=alive)
    if (.not. alive) then
        write(*,*) "Control file: """//trim(ConFile)//""" doesn't exist!"
        pause
        stop
    end if
    ConID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=ConID, file=ConFile)
    
    do while(.true.)
        read(unit=ConID,fmt='(A)',end=200) line
        if (index(line,"!") /= 0)  line=line(1:index(line,"!")-1)
        if (len_trim(line)==0) cycle
        if (index(line,"day") /= 0)   then
            read(line,*) temp, int_year, int_doy
        elseif (index(line,"eopfile") /= 0)   then
            read(line,"(12X,A)") EOPFile
        elseif (index(line,"antfile") /= 0)   then
            read(line,"(12X,A)") AntFile
        elseif (index(line,"planetfile") /= 0)   then
            read(line,"(12X,A)") PlanetFile
        elseif (index(line,"tidefile") /= 0)   then
            read(line,"(12X,A)") TideFile
        elseif (index(line,"glofrefile") /= 0)   then
            read(line,"(12X,A)") GloFreFile
        elseif (index(line,"gpt2griddir") /= 0)   then
            read(line,"(12X,A)") GPT2GridDir
        elseif (index(line,"dcbfile") /= 0)   then
            read(line,"(12X,A)") DCBFile
        elseif (index(line,"p1c1file") /= 0)   then
            read(line,"(12X,A)") P1C1File
        elseif (index(line,"obstype") /= 0)   then
            read(line,*) temp, ObsType
        elseif (index(line,"IORQ") /= 0)   then
            read(line,*) temp, IorQ
        elseif (index(line,"obsdir") /= 0)   then
            read(line,"(12X,A)") ObsDir
        elseif (index(line,"orbit") /= 0)   then
            read(line,*) temp, Orbit
        elseif (index(line,"navdir") /= 0)   then
            read(line,"(12X,A)") NavDir
        elseif (index(line,"sp3dir") /= 0)   then
            read(line,"(12X,A)") sp3dir
        elseif (index(line,"clock") /= 0)   then
            read(line,*) temp, Clk
        elseif (index(line,"clkdir") /= 0)   then
            read(line,"(12X,A)") clkdir
        elseif (index(line,"coortable") /= 0)   then
            read(line,"(12X,A)") CoorTable
        elseif (index(line,"interval") /= 0)   then
            read(line,*) temp, Interval
        elseif (index(line,"clktype") /= 0)   then
            read(line,*) temp, clktype
        elseif (index(line,"limele") /= 0)   then
            read(line,*) temp, LimEle
        elseif (index(line,"limsnr") /= 0)   then
            read(line,*) temp, LimSNR
        elseif (index(line,"maxpdop") /= 0)   then
            read(line,*) temp, MaxPDOP
        elseif (index(line,"tropsphere") /= 0)   then
            read(line,*) temp, TropLen
        elseif (index(line,"trop_model") /= 0)   then
            read(line,*) temp, cdattype, cztd, cmap
        elseif (index(line,"pos_state") /= 0)   then
            read(line,*) temp, str_pos_state
             if ( (index(str_pos_state,"K") /=0) .or. (index(str_pos_state,"k") /=0)  )then
                Pos_State="K"
            elseif ( (index(str_pos_state,"F") /=0) .or. (index(str_pos_state,"f") /=0)  )then
                Pos_State="F"
            end if
        elseif (index(line,"start_time") /= 0)   then
            read(line,*) temp, year,mon,day, hour, min, sec
            call UTC2GPST(year, mon, day,  hour, min, sec, GPSweek_st, GPSsec_st)
        elseif (index(line,"end_time") /= 0)   then
            read(line,*) temp, year,mon,day, hour, min, sec
            call UTC2GPST(year, mon, day,  hour, min, sec, GPSweek_end, GPSsec_end)
        elseif (index(line,"diff_mod") /= 0)   then
            read(line,*) temp, diffmod
        elseif (index(line,"freq_comb") /= 0)   then
            read(line,*) temp, freq_comb
        elseif (index(line,"cycleslip") /= 0)   then
            read(line,*) temp, csmethod
        elseif (index(line,"obscombine") /= 0)   then
            read(line,*) temp, ObsCombine
             if  ((index(ObsCombine,"PC") ==0) .and. (index(ObsCombine,"G1") ==0) .and. (index(ObsCombine,"G2") ==0)  .and. (index(ObsCombine,"G3") ==0) &
                  .and. (index(ObsCombine,"P1") ==0) .and. (index(ObsCombine,"P2") ==0) .and. (index(ObsCombine,"P3") ==0))  then
                write(*,*) 'Obs combination incorrect. Please check.'
                pause 
                stop
            end if
        elseif (index(line,"proc_mod") /= 0)   then
            read(line,*) temp, proc_mod
        elseif (index(line,"dd_coe") /= 0)   then
            read(line,*) temp, a1,a2,b1,b2
        elseif (index(line,"partial_ar") /= 0)   then
            read(line,*) temp, str_partial_ar
            if ( (index(str_partial_ar,"N") /=0) .or.  (index(str_partial_ar,"n") /=0) )then
                partial_ar=.false.
            else
                read(line,*) temp, str_partial_ar, parARnum
            end if
        elseif (index(line,"ratio") /= 0)   then
            read(line,*) temp, minratio
        elseif (index(line,"ar_mode") /= 0)   then
            read(line,*) temp, ar_mode
        elseif (index(line,"fixholdele") /= 0)   then
            read(line,*) temp, FixEle, HoldEle
        elseif (index(line,"est_iono") /= 0)   then
            read(line,*) temp, temp
            if ( (index(temp,"Y") /=0) .or. (index(temp,"y") /=0)  )then
                If_Est_Iono=.true.
            end if
        elseif (index(line,"est_wl_amb") /= 0)   then
            read(line,*) temp, temp
            if ( (index(temp,"Y") /=0) .or. (index(temp,"y") /=0)  )then
                If_Est_WL=.true.
            end if
        elseif (index(line,"ionocompensate") /= 0)   then
            read(line,*) temp, str_IonoCompensate
            if ( (index(str_IonoCompensate,"Y") /=0) .or.  (index(str_IonoCompensate,"y") /=0) )then
                If_IonoCompensate=.true.
            end if
        elseif (index(line,"pcorfile") /= 0)   then
            read(line,"(12X,A)")  pcorfile
        elseif (index(line,"orbcorrfile") /= 0)   then
            read(line,"(12X,A)")  orbcorrfile
        elseif (index(line,"zonecorrfile") /= 0)   then
            read(line,"(12X,A)")  zonecorrfile
        elseif (index(line,"bebfile") /= 0)   then
            read(line,"(12X,A)")  BEBFile
        elseif (index(line,"iondir") /= 0)   then
            read(line,"(12X,A)")   IonDir
        elseif (index(line,"iontype") /= 0)   then
            read(line,*) temp, iontype
        elseif (index(line,"delay") /= 0)   then
            read(line,*) temp, delay0
            delay=delay*60.d0  ! minute to second
        elseif (index(line,"smooth") /= 0)   then
            read(line,*) temp, Var_smooth, Smooth_Method, Smooth_Combine, Smooth_Time
            if (Smooth_Time>Interval*59) Smooth_Time=Interval*59
       elseif (index(line,"combination") /= 0)   then
            read(line,*) temp, str_Combination,sigPC
            Combination(1)=.true.
            if (index(str_Combination,"PC") ==0) then
                write(*,*) "PC combination not found. Please check."
                pause
                stop
            end if
            if (index(str_Combination,"LC") /=0) then
                Combination(2)=.true.
                read(line,*) temp, str_Combination,sigPC, sigLC
            end if
            if (index(str_Combination,"DP") /=0) then
                Combination(3)=.true.
                read(line,*) temp, str_Combination,sigPC, sigLC, sigDP
            end if
        elseif (index(line,"systemused") /= 0)   then
            read(line,"(12X,A20)") str_SystemUsed
            SystemUsed(1)=.false.
            if ((index(str_SystemUsed,"G") /=0) .or. (index(str_SystemUsed,"g") /=0)  )then
                SystemUsed(1)=.true.
            end if
            if ((index(str_SystemUsed,"R") /=0) .or. (index(str_SystemUsed,"r") /=0)  )then
                SystemUsed(2)=.true.
            end if
            if ((index(str_SystemUsed,"C") /=0) .or. (index(str_SystemUsed,"c") /=0)  )then
                SystemUsed(3)=.true.
            end if
            if ((index(str_SystemUsed,"E") /=0) .or. (index(str_SystemUsed,"e") /=0)  )then
                SystemUsed(4)=.true.
            end if
            if ((index(str_SystemUsed,"J") /=0) .or. (index(str_SystemUsed,"j") /=0)  )then
                SystemUsed(5)=.true.
            end if
            if ((index(str_SystemUsed,"I") /=0) .or. (index(str_SystemUsed,"i") /=0)  )then
                SystemUsed(6)=.true.
            end if
        elseif (index(line,"admethod") /= 0)   then
            read(line,*) temp, admethod
            if ((index(admethod,"LS") ==0) .and. (index(admethod,"KF") ==0)  )then
                write(*,*) "Adjustment method incorrect. Please check."
                pause
                stop
            end if
        elseif (index(line,"isb") /= 0)   then
            read(line,*) temp, str_ISB
            if ( (index(str_ISB,"NO") /=0) .or.  (index(str_ISB,"no") /=0) )then
                If_ISB=.false.
            else 
                read(line,*) temp, str_ISB, isb_mode
                if ( (index(isb_mode,"RW") /=0) .or. (index(isb_mode,"rw") /=0) ) then
                    read(line,*) temp, str_ISB,  isb_mode, isb_step
                elseif ( (index(isb_mode,"CV") /=0) .or. (index(isb_mode,"cv") /=0) ) then
                    
                elseif ( (index(isb_mode,"WN") /=0) .or. (index(isb_mode,"wn") /=0) ) then
                    
                else
                    write(*,*) 'Can not identify ISB mode:  ', str_ISB
                    pause
                    stop
                end if
            end if
        elseif (index(line,"ifb") /= 0)   then
            read(line,*) temp, str_IFB
            if ( (index(str_IFB,"YES") /=0) .or.  (index(str_IFB,"yes") /=0) )then
                If_IFB=.true.
                read(line,*) temp, str_IFB, IFB_Mode
                if ( (index(IFB_Mode,"FD") ==0) .and. (index(IFB_Mode,"fd") ==0) .and. (index(IFB_Mode,"LM") ==0) .and. (index(IFB_Mode,"lm") ==0) ) then
                    write(*,*) 'Can not identify IFB mode:  ', IFB_Mode
                    pause
                    stop
                end if
            end if
        elseif (index(line,"outdir") /= 0)   then
            read(line,"(12X,A)") OutDir
        end if
    end do
    200 close(ConID)
!    idx=len_trim(sp3dir)    ! If the gbm orbit, use the gbm PCO since DOY 2014 197
!    if ( (index(sp3dir(idx-2:idx),'gbm') /=0) .and. ((int_year>=2015) .or. ((int_year==2014) .and. (int_doy>197))) ) then
!        AntFile=trim(AntFile)//'gbm'
!    elseif (index(sp3dir(idx-2:idx),'wum') /=0) then ! If the wum orbit, use the wum PCO
!        AntFile=trim(AntFile)//'wum'
!    end if
    if (FixEle<LimEle) FixEle=LimEle
    
    if (SystemUsed(1)) then
        if (Combination(1)) SatNum=SatNum+GNum
        ParaNum=ParaNum+1
        INT_SystemUsed(1)=1
    else
        GNum=0
    end if
    if (SystemUsed(2)) then
        if (Combination(2)) SatNum=SatNum+RNum
        ParaNum=ParaNum+1
        INT_SystemUsed(2)=1
        if (If_IFB) then
            GloParaNum=14
            if ( (index(IFB_Mode,"FD") /=0) .or. (index(IFB_Mode,"fd") /=0) ) then
                IFB_Mode='FD'
                ParaNum=ParaNum+14  ! Frequency depend
            elseif ( (index(IFB_Mode,"LM") /=0) .or. (index(IFB_Mode,"lm") /=0) ) then
                IFB_Mode='LM'
                ParaNum=ParaNum+1  ! Linear model
            end if
        end if
    else
        RNum=0
    end if
    if (SystemUsed(3)) then
        if (Combination(2)) SatNum=SatNum+CNum
        ParaNum=ParaNum+1
        INT_SystemUsed(3)=1
    else
        CNum=0
    end if
    if (SystemUsed(4)) then       ! GALILEO
        ParaNum=ParaNum+1
        INT_SystemUsed(4)=1
        if (Combination(2)) SatNum=SatNum+NumE
    else
        NumE=0
    end if
    if (SystemUsed(5)) then       ! QZSS
        ParaNum=ParaNum+1
        INT_SystemUsed(5)=1
        if (Combination(2)) SatNum=SatNum+JNum
    else
        JNum=0
    end if
    if (SystemUsed(6)) then       ! IRNSS
        ParaNum=ParaNum+1
        INT_SystemUsed(6)=1
        if (Combination(2)) SatNum=SatNum+INum
    else
        INum=0
    end if
    if (If_Est_Iono .and. (index(ObsCombine,"P1") /=0 .or. index(ObsCombine,"P2") /=0 .or. index(ObsCombine,"P3") /=0)) then       ! If estimate ionosphere
        ParaNum=ParaNum+SatNum
    end if

    if (GPSweek_st/=0) then
        n=((GPSweek_end-GPSweek_st)*604800.d0+GPSsec_end-GPSsec_st)/7200+4
        if (n<0) then
            write(*,*) "输入的起始时间大于终止时间"
            pause
            stop
        elseif (n<24) then
            n=24
        end if
    else
        n=25  !如果没有起止时间
        call UTC2GPST(int_year, 1, int_doy,  0, 0, 0.d0, GPSweek_st, GPSsec_st)
        call UTC2GPST(int_year, 1, int_doy,  23, 59, 59.d0, GPSweek_end, GPSsec_end)
    end if
    do PRN=1,GNum  ! GPS ephemeris
        allocate(NavData(PRN)%Nav(n))
    end do
    do PRN=1,RNum  ! GLONASS ephemeris
        allocate(NavData_R(PRN)%Nav(n*4+2))
    end do
    do PRN=1,CNum  ! BD ephemeris
        allocate(NavData_C(PRN)%Nav(2*n))
    end do
    do PRN=1,NumE  ! GALILEO ephemeris
        allocate(NavData_E(PRN)%Nav(24*n))
    end do
    do PRN=1,JNum  ! QZSS ephemeris
        allocate(NavData_J(PRN)%Nav(2*n))
    end do
    do PRN=1,INum  ! IRNSS ephemeris
        allocate(NavData_I(PRN)%Nav(8*n))
    end do
    return
end subroutine