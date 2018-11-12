! ============== Get_DCB =================
!
! PURPOSE:
!        1.  Get P1C1 DCBs of GPS and GLONASS from DCB files
!        2. Get BDS DCB from MGEX DCB file or BDS DCBIQ.dat
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================  End of header  =================

subroutine Get_DCB
use MOD_VAR
use MOD_FileID
use MOD_FileDir
implicit none
    logical   :: alive, aliveP1C1
    character(200) :: line
    integer   :: error,PRN
    character(1) :: System
    real(8)   :: temp(5), val
    integer :: year_st, year_end, doy_st, doy_end


    ! Read P1P2 DCB File
    DCB=0.d0
    inquire(file=P1C1File,exist=aliveP1C1)
    if (.not. aliveP1C1) then
        write(*,*) "P1C1 file: """//trim(P1C1File)//""" doesn't exist. It will set as zero or use Multi-GNSS DCB(*.bsx)."
        !pause
        goto 100
    end if
    P1C1ID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=P1C1ID, file=P1C1File)
    
    do while(.true.)
        read(unit=P1C1ID,fmt='(A)', iostat=error) line
        if (error /=0) exit
        if (index(line,"***") /= 0)  exit
    end do 

    do while(.true.)
        read(unit=P1C1ID,fmt='(A)', iostat=error) line
        if (error /=0) exit
        if (len_trim(line)==0) cycle
        read(line(1:1),*) System
        read(line(2:200),*) PRN,temp(1:2)
        if (System=='R') then
            if (.not. SystemUsed(2)) cycle
            PRN=PRN+GNum
        elseif  (System=='G') then
            if (.not. SystemUsed(1)) cycle
        end if
        DCB(PRN)=temp(1)*1d-9
    end do

    close(P1C1ID)

    ! Read DCB File
    100 DCBBSX=0.d0
    inquire(file=DCBFile,exist=alive)
    if (.not. alive) then
        if ( (index(Clk,"CLK")/=0) .and. ( (index(ObsCombine,"PC")==0) .or. (index(freq_comb,'L1L2')==0) ) ) then  ! If precise clock and not L1L2 combination
            write(*,*) "No Multi-GNSS DCB exist. It will affect the real clock."
            pause
        elseif ( (index(Clk,"BRD")/=0) .and. (SystemUsed(2)) .and. ( (index(ObsCombine,"P1")==0) .and. (index(ObsCombine,'G1')==0) ) ) then  ! If precise clock and not L1L2 combination
            write(*,*) "No Multi-GNSS DCB exist. It will affect the real clock of GLONASS."
!            pause
        end if
        return
    end if
    DCBID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=DCBID, file=DCBFile)
    
    if (index(DCBFile,'DCBIQ')/=0) then  ! BDS DCB file
        do while(.true.)
            read(unit=DCBID,fmt='(A)', iostat=error) line
            if (error /=0) exit
            read(line,*) PRN,temp(1:5)
            if (PRN>35) cycle
            temp=temp*1d-9
            DCBIQ(1,PRN)=temp(2)
            DCBIQ(2,PRN)=temp(2)-temp(1)
            DCBIQ(3,PRN)=temp(4)
            DCBIQ(4,PRN)=temp(4)-temp(3)
            DCBIQ(5,PRN)=temp(5)
        end do
    else                         ! MGEX DCB
        do while(.true.)
            read(unit=DCBID,fmt='(A)', iostat=error) line
            if (error /=0) exit
            if (index(line, "*BIAS") /= 0)  exit
        end do 

        do while(.true.)
            read(unit=DCBID,fmt="(A)") line
            if ((index(line,"-BIAS/SOLUTION") /= 0) .or. (index(line(16:19),"    ")==0) ) then
                exit
            elseif (index(line,"OSB") /= 0) then
                backspace(DCBID)
                return
            elseif ((index(line(1:4),"DSB") /= 0) .or. (index(line(1:4),"DCB") /= 0)) then
                if (index(line(1:4),"DSB") /= 0) then
                    read(line(38:39), fmt="(I2)") year_st
                    read(line(41:43), fmt="(I3)") doy_st
                    read(line(53:54), fmt="(I2)") year_end
                    read(line(56:58), fmt="(I3)") doy_end
                elseif (index(line(1:4),"DCB") /= 0) then
                    read(line(41:42), fmt="(I2)") year_st
                    read(line(44:46), fmt="(I3)") doy_st
                    read(line(54:55), fmt="(I2)") year_end
                    read(line(57:59), fmt="(I3)") doy_end
                end if
                if ( ((int_year-2000-year_st)*365+int_doy-doy_st<300) .and. ((int_year-2000-year_st)*365+int_doy-doy_st>=0) ) then
                    read(line(13:14), fmt="(I2)") PRN
                    read(line(12:12), fmt="(A1)") System
                    read(line(81:92), fmt="(F12.4)") val
                    if (System=='G') then
                        if (.not. SystemUsed(1)) cycle
                        if (index(line,"C1C  C2W") /= 0) then  ! C1P2
                            DCBBSX(1,PRN)=val*1e-9
                        elseif (index(line,"C1C  C5Q") /= 0) then   ! C1P5
                            DCBBSX(2,PRN)=val*1e-9
                        elseif (index(line,"C1C  C1W") /= 0) then  ! C1P1
                            DCBBSX(3,PRN)=val*1e-9
                            if (.not. aliveP1C1) then
                                DCB(PRN)= -val*1e-9  
                                ! IGG/DLR DCB is oppisite of CODE and ISCL1CA(CNAV) 
                                ! Reference: 王宁波, 袁运斌, 张宝成,等. GPS民用广播星历中ISC参数精度分析及其对导航定位的影响[J]. 测绘学报, 2016, 45(8):919-928.
                                ! Steigenberger P, Montenbruck O, Hessel U. Performance Evaluation of the Early CNAV Navigation Message[J]. Navigation, 2015, 62(3):219-228
                            end if
                        else
                        end if
                    elseif (System=='R') then
                        if (.not. SystemUsed(2)) cycle
                        if (index(line,"C1C  C2P") /= 0) then  ! C1P2
                            DCBBSX(1,PRN+GNum)=val*1e-9
                        elseif (index(line,"C1P  C2P") /= 0) then  ! P1P2
                            DCBBSX(2,PRN+GNum)=val*1e-9
                        elseif (index(line,"C1C  C1P") /= 0) then  ! C1P1
                            DCBBSX(3,PRN+GNum)=val*1e-9
                            DCB(PRN+GNum)= -val*1e-9  
                        elseif (index(line,"C1C  C2C") /= 0) then  ! C1C2C
                            DCBBSX(4,PRN+GNum)=val*1e-9
                        end if
                    elseif (System=='E') then
                        if (.not. SystemUsed(4)) cycle
                        if (index(line,"C1C  C5Q") /= 0) then  ! C1P2
                            DCBBSX(1,PRN+GNum+RNum+CNum)=val*1e-9
                        elseif (index(line,"C1C  C7Q") /= 0) then   ! C1P3
                            DCBBSX(2,PRN+GNum+RNum+CNum)=val*1e-9
                        end if
                    elseif (System=='J') then
                        if (.not. SystemUsed(4)) cycle
                        if (index(line,"C1C  C2L") /= 0) then  ! C1P2
                            DCBBSX(1,PRN+GNum+RNum+CNum+NumE)=val*1e-9
                        elseif (index(line,"C1C  C5Q") /= 0) then   ! C1P5
                            DCBBSX(2,PRN+GNum+RNum+CNum+NumE)=val*1e-9
                        elseif (index(line,"C1C  C1X") /= 0) then   ! C1P1
                            DCBBSX(3,PRN+GNum+RNum+CNum+NumE)=val*1e-9
                            if (DCB(PRN+GNum+RNum+CNum+NumE)==0.d0) then
                                DCB(PRN+GNum+RNum+CNum+NumE)= -val*1e-9  
                            end if
                        elseif (index(line,"C1X  C2X") /= 0) then   ! C1X C2X
                            DCBBSX(4,PRN+GNum+RNum+CNum+NumE)=val*1e-9
                        end if
                    elseif (System=='C') then
                        if (.not. SystemUsed(3)) cycle
                        if (index(line,"C2I  C7I") /= 0) then  ! B1B2
                            DCBIQ(5,PRN)=val*1e-9
                            DCBBSX(1,PRN+GNum+RNum)=val*1e-9
                        elseif (index(line,"C2I  C6I") /= 0) then   ! B1B3
                            DCBIQ(1,PRN)=val*1e-9
                            DCBIQ(3,PRN)=val*1e-9
                            DCBBSX(2,PRN+GNum+RNum)=val*1e-9
                        elseif (index(line,"C7I  C6I") /= 0) then  ! B2B3
                            DCBIQ(2,PRN)=val*1e-9
                            DCBIQ(4,PRN)=val*1e-9
                        end if
                    end if
                end if ! if ( ((int_year-20
            end if
        end do
        
        if (all(DCBBSX==0.d0)) then
            write(*,*) "Date of Multi-GNSS DCB difference exceeds 100 days, it will set as 0."
!            pause
        end if

        ! Change GPS DCB(C1C,C2W) to DCB(C1W,C2W)
        ! Change GLONASS DCB(C1C,C2P) to DCB(C1P,C2P)
        do PRN=1,GNum+RNum
            DCBBSX(1,PRN)=DCBBSX(1,PRN)-DCBBSX(3,PRN)
        end do

        ! If only B1B2 and B1B3
        do PRN=1,15
            if ((DCBIQ(5,PRN)/=0.d0) .and. (DCBIQ(3,PRN)/=0.d0) .and. (DCBIQ(4,PRN)==0.d0)) then
                DCBIQ(4,PRN)=DCBIQ(3,PRN)-DCBIQ(5,PRN)
                DCBIQ(2,PRN)=DCBIQ(4,PRN)
                DCBIQ(5,PRN)=0.d0
            end if
        end do

    end if 

    close(DCBID)

    return

end subroutine

