!============ Get_CNESBias ==========
!
! PURPOSE:
!      This function is used to get CNES code 
!   and phase bias for PPP AR
!
! INPUTS:
!       ObsWeek           
!       ObsSow
!
! WRITTEN BY: Yize Zhang, zhyize@163.com
!===========End of Header==========

subroutine Get_CNESBias(ObsWeek, ObsSow) 
use MOD_Bias
use MOD_FileID
use MOD_VAR
implicit none
    integer :: ObsWeek
    real(8) :: ObsSow
    ! Local variables
    integer :: error
    character(100) :: line
    integer :: year_st, year_end, doy_st, doy_end, sod_st, sod_end
    integer :: Week_St, Week_End
    real(8) :: Sow_St,Sow_End
    integer :: PRN
    character(1) :: System
    real(8) :: val

    do while(.true.)
        read(DCBID,fmt='(A)', iostat=error) line
        if (index(line,"-BIAS/SOLUTION") /= 0) then
            exit
        end if
        
        read(line(36:39), fmt="(I4)") year_st
        read(line(41:43), fmt="(I3)") doy_st
        read(line(45:49), fmt="(I5)") sod_st
        read(line(51:54), fmt="(I4)") year_end
        read(line(56:58), fmt="(I3)") doy_end
        read(line(60:64), fmt="(I5)") sod_end

        call UTC2GPST(year_st, 1, doy_st,  0, 0, sod_st+0.001d0, Week_St, Sow_St)
        call UTC2GPST(year_end, 1, doy_end,  0, 0, sod_end+0.001d0, Week_End, Sow_End)

        if ( (Week_St-ObsWeek)*604800.d0+(Sow_St-ObsSow)>0.01d0) then
            backspace(DCBID)    ! Íí£¬ÍË³ö
            exit
        elseif ( (Week_End-ObsWeek)*604800.d0+(Sow_End-ObsSow)<-1800.d0) then
            cycle ! Ôç£¬¼ÌÐø
        end if
        
        Bias%week=Week_St
        Bias%sow=Sow_St-0.001d0
        read(line(13:14), fmt="(I2)") PRN
        read(line(12:12), fmt="(A1)") System
        read(line(81:91), fmt="(F11.4)") val
        val=val*1.d-9  ! Change to second

        if (System=='G') then
            if (.not. SystemUsed(1)) cycle
            if (index(line,"C1W") /= 0) then
                Bias%P1(PRN)=val
            elseif (index(line,"C2W") /= 0) then 
                Bias%P2(PRN)=val
            elseif (index(line,"C5Q") /= 0) then 
                Bias%P3(PRN)=val
            elseif (index(line,"L1C") /= 0) then 
                Bias%L1(PRN)=val
            elseif (index(line,"L2W") /= 0) then 
                Bias%L2(PRN)=val
            elseif (index(line,"L5I") /= 0) then
                Bias%L3(PRN)=val
            end if
        elseif (System=='R') then
            if (.not. SystemUsed(2)) cycle
            PRN=PRN+GNum
            if (index(line,"C1P") /= 0) then
                Bias%P1(PRN)=val
            elseif (index(line,"C2P") /= 0) then
                Bias%P2(PRN)=val
            end if
        elseif (System=='C') then
            if (.not. SystemUsed(3)) cycle
            PRN=PRN+GNum+RNum
            if (index(line,"C1I") /= 0 .or. index(line,"C2I") /= 0) then
                Bias%P1(PRN)=val
            elseif (index(line,"C7I") /= 0) then
                Bias%P2(PRN)=val
            elseif (index(line,"C6I") /= 0) then
                Bias%P3(PRN)=val
            elseif (index(line,"L1I") /= 0 .or. index(line,"L2I") /= 0) then
                Bias%L1(PRN)=val
            elseif (index(line,"L7I") /= 0) then
                Bias%L2(PRN)=val
            elseif (index(line,"L6I") /= 0) then
                Bias%L3(PRN)=val
            end if
        elseif (System=='E') then
            if (.not. SystemUsed(4)) cycle
            PRN=PRN+GNum+RNum+CNum
            if (index(line,"C1X") /= 0) then
                Bias%P1(PRN)=val
            elseif (index(line,"C5X") /= 0) then
                Bias%P2(PRN)=val
            elseif (index(line,"C7X") /= 0) then
                Bias%P3(PRN)=val
            elseif (index(line,"L1X") /= 0) then
                Bias%L1(PRN)=val
            elseif (index(line,"L5X") /= 0) then
                Bias%L2(PRN)=val
            elseif (index(line,"L7X") /= 0) then
                Bias%L3(PRN)=val
            end if
        end if
    end do
    
    return
end subroutine