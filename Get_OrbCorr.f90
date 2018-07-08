!==================== Get_OrbCorr==================
!
! PURPOSE:
!      This function is used to get Orbit Correction.
!      Adapted from Get_ESC.
!
! INPUTS:
!       ObsWeek        week     
!       ObsSow         second of week
!       IQ             0： I band;
!                      3:  Q band
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!==================== End of Header ==============

subroutine Get_OrbCorr(ObsWeek, ObsSow, IQ) 
use MOD_ESC
use MOD_FileID
implicit none
    integer :: ObsWeek, Week
    real(8) :: ObsSow, Sow, dt
    integer(1) :: PRN, error, IQ, IQt
    character(100) :: line, temp

    OrbCorr=0.d0
    dt=0.d0
    if (IQ==0) dt=1.d0
    Week=ObsWeek;Sow=ObsSow
    do while(.true.)
        read(OrbCorrID,fmt='(a)', iostat=error)  line
        if (error /=0) exit
        if (index(line,"*****") /= 0)   cycle
        if (index(line,"week0") /= 0)   cycle
        read(line, *) week, sow, PRN , IQt
        week=week+1356
        if ( (Week-ObsWeek)*604800.d0+(Sow+dt-ObsSow)>0.d0) then
            backspace(OrbCorrID)   ! 晚，退出
            exit
        end if
        if ( (Week-ObsWeek)*604800.d0+(Sow+dt-ObsSow)<0.d0) then
            cycle  ! 早，继续
        end if
        if (IQt==IQ) then
            read(line(24:32),*) OrbCorr(PRN+56,1)
            read(line(54:63),*) OrbCorr(PRN+56,2)
            read(line(84:93),*) OrbCorr(PRN+56,3)
            OrbCorrSow=sow + dt
        end if

    end do

    if ( (abs(Week-ObsWeek)*604800.d0+Sow-ObsSow)> 43200.d0) then
        write(*,*) 'X37a文件中的时间与观测时间差异超过12小时，请检查。'
        pause
    end if

end subroutine