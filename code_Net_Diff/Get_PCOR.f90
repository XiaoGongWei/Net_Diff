!==================== Get_PCOR==================
!
! PURPOSE:
!      This function is used to get PCOR.
!      Adapted from Get_ESC.
!
! INPUTS:
!       ObsWeek           
!       ObsSow
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!==================== End of Header ==============

subroutine Get_PCOR(ObsWeek, ObsSow, IQ) 
use MOD_ESC
use MOD_FileID
use MOD_VAR
implicit none
    integer :: ObsWeek, Week
    real(8) :: ObsSow, Sow, dt
    integer(1) :: PRN, error, IQ, IQt
    character(100) :: line, temp
    real(8) :: pcor1, pcor2

    ESC=0.d0
    dt=0.d0
    if (IQ==0) dt=1.d0
    Week=ObsWeek;Sow=ObsSow
    do while(.true.)
        read(PcorID,fmt='(a)', iostat=error)  line
        if (error /=0) exit
        if (index(line,"*****") /= 0)   cycle
        read(line, *) week, sow, PRN , IQt
        week=week+1356
        if ( (Week-ObsWeek)*604800.d0+(Sow+dt-ObsSow)>0.d0) then
            backspace(PcorID)   ! 晚，退出
            exit
        end if
        if ( (Week-ObsWeek)*604800.d0+(Sow+dt-ObsSow)<0.d0) then
            cycle  ! 早，继续
        end if
        if (IQt==IQ) then
            if ( (freq_comb=="L1L2") .or. (proc_mod==3) ) then
                read(line(26:36),*) pcor1
                read(line(50:60),*) pcor2
                ESC(PRN+GNum+RNum)=(152.6d0**2*pcor1-118.d0**2*pcor2)/(152.6d0+118.d0)/(152.6d0-118.d0)
            elseif (freq_comb=="L1L3") then
                read(line(26:36),*) pcor1
                read(line(74:84),*) pcor2
                ESC(PRN+GNum+RNum)=(152.6d0**2*pcor1-124.d0**2*pcor2)/(152.6d0+124.d0)/(152.6d0-124.d0)
            elseif (freq_comb=="L2L3") then
                read(line(50:60),*) pcor1
                read(line(74:84),*) pcor2
                ESC(PRN+GNum+RNum)=(118.0d0**2*pcor1-124.d0**2*pcor2)/(118.0d0+124.d0)/(118.0d0-124.d0)
            end if
            PCORSow=sow + dt
        end if
        ESC(PRN+GNum+RNum)=-ESC(PRN+GNum+RNum)
    end do
!    ESC=0.d0

    if ( (abs(Week-ObsWeek)*604800.d0+Sow-ObsSow)> 43200.d0) then
        write(*,*) 'PCOR文件中的时间与观测时间差异超过12小时，请检查。'
        pause
    end if

end subroutine