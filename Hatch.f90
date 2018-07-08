! ==================  Hatch_Filter  ================
!
! PURPOSE:
!            Smooth the pseudo-range combined with phase observations
!            by using Hatch flitering method.
!
! INPUTS:
!          PRN               PRN of the satellite
!          P1,P2             Pseudorange observation in L1 & L2 band, in meters
!          Range            Pseudo-range before smoothing
!          L1,L2              Phase observation in L1 & L2 band, in meters
!          Ele                 Elevation of the satellite
!          k                     station order
!          Slip                 Flag of the cycle slip of phase observations
! 
! OUTPUT:
!          Range            Smoothed pseudo-range
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ===============  End of Header  ==================

subroutine Hatch_Filter(PRN, P1, P2, Range, L1, L2, Ele, k, Slip)
use MOD_FileID
use MOD_Var
use MOD_STA
implicit none
    ! Intent
    integer     :: PRN
    integer(1) :: Slip
    real(8) :: P1,P2, L1,L2, Ele
    integer :: k
    ! Local variables
    real(8) :: Range
    real(8) :: Phase
    real(8) :: Mp   ! Multipath error
    integer :: i


    if (Slip==1) then  ! If cycle slip, initialization
        STA%STA(k)%Pre(PRN)%n=0.d0
        STA%STA(k)%Pre(PRN)%Mp1=0.d0
        STA%STA(k)%Pre(PRN)%Mp2=0.d0
        STA%STA(k)%Pre(PRN)%HisRange=99.d0
        STA%STA(k)%Pre(PRN)%HisPhase=99.d0
    end if
    
    Phase=0.d0
    if ( (L1/=0.d0) .and. (L2/=0.d0)) then
        Phase=(f1*f1*L1-f2*f2*L2)/(f1+f2)/(f1-f2)
     end if
    
    ! Hatch滤波进行相位平滑伪距
    if (index(Smooth_Combine,"PC")/=0) then
        if ((STA%STA(k)%Pre(PRN)%Range/=0.d0) .and. (STA%STA(k)%Pre(PRN)%PrePhase/=0.d0) .and.(Phase/=0.d0) .and.(Range/=0.d0) ) then
            STA%STA(k)%Pre(PRN)%n=STA%STA(k)%Pre(PRN)%n+1.d0 ! dsind(Ele)
            Range=Range*1.d0/STA%STA(k)%Pre(PRN)%n+(1.d0-1.d0/STA%STA(k)%Pre(PRN)%n)*(STA%STA(k)%Pre(PRN)%Range+(Phase-STA%STA(k)%Pre(PRN)%PrePhase) )   ! 平滑后的无电离层组合
        end if
    ! CNMC进行相位平滑伪距
    else if (index(Smooth_Combine,"P1P2")/=0) then
        if ( (L1/=0.d0) .and. (L2/=0.d0)  .and. (P1/=0.d0) .and. (P2/=0.d0) ) then    ! CNMC +STA%STA(k)%Pre(PRN)%Mp1
            STA%STA(k)%Pre(PRN)%n=STA%STA(k)%Pre(PRN)%n+1.d0 ! dsind(Ele)
            STA%STA(k)%Pre(PRN)%CBias1=(P1-L1 -2*f2*f2/(f1+f2)/(f1-f2)*(L1-L2))*1.d0/STA%STA(k)%Pre(PRN)%n+(1.d0-1.d0/STA%STA(k)%Pre(PRN)%n)*(STA%STA(k)%Pre(PRN)%CBias1)
            STA%STA(k)%Pre(PRN)%Mp1=L1+2*f2*f2/(f1+f2)/(f1-f2)*(L1-L2)+STA%STA(k)%Pre(PRN)%CBias1-P1
            P1=STA%STA(k)%Pre(PRN)%Mp1+P1
            STA%STA(k)%Pre(PRN)%CBias2=(P2-L2-2*f1*f1/(f1+f2)/(f1-f2)*(L1-L2))*1.d0/STA%STA(k)%Pre(PRN)%n+(1.d0-1.d0/STA%STA(k)%Pre(PRN)%n)*(STA%STA(k)%Pre(PRN)%CBias2)
            STA%STA(k)%Pre(PRN)%Mp2=L2+2*f1*f1/(f1+f2)/(f1-f2)*(L1-L2)+STA%STA(k)%Pre(PRN)%CBias2-P2
            P2=STA%STA(k)%Pre(PRN)%Mp2+P2
        end if
        if ( (P1/=0.d0) .and. (P2/=0.d0) ) then
            Range=(f1*f1*P1-f2*f2*P2)/(f1+f2)/(f1-f2) ! Ionospheric-free combination pseudo-range after smoothing
        end if
    elseif ( (index(Smooth_Combine,"P1")/=0) .or. (index(Smooth_Combine,"P2")/=0) .or. (index(Smooth_Combine,"P3")/=0) ) then
        if (index(Smooth_Combine,"P1")/=0) then
            STA%STA(k)%Pre(PRN)%HisRange=eoshift(STA%STA(k)%Pre(PRN)%HisRange,shift=1,boundary=99.d0)
            STA%STA(k)%Pre(PRN)%HisPhase=eoshift(STA%STA(k)%Pre(PRN)%HisPhase,shift=1,boundary=99.d0)
            if ((P1/=0.d0) .and. (L1/=0.d0)) then
                STA%STA(k)%Pre(PRN)%HisRange(60)=P1
                STA%STA(k)%Pre(PRN)%HisPhase(60)=L1
                Range=P1
                Phase=L1
            end if
        elseif ( (index(Smooth_Combine,"P2")/=0)  .and. (freq_comb=='L2L3')) then
            STA%STA(k)%Pre(PRN)%HisRange=eoshift(STA%STA(k)%Pre(PRN)%HisRange,shift=1,boundary=99.d0)
            STA%STA(k)%Pre(PRN)%HisPhase=eoshift(STA%STA(k)%Pre(PRN)%HisPhase,shift=1,boundary=9.d0)
             if ((P1/=0.d0) .and. (L1/=0.d0)) then
                STA%STA(k)%Pre(PRN)%HisRange(60)=P1
                STA%STA(k)%Pre(PRN)%HisPhase(60)=L1
                Range=P1
                Phase=L1
            end if
        elseif ( (index(Smooth_Combine,"P2")/=0)  .or. (index(Smooth_Combine,"P3")/=0) ) then
            STA%STA(k)%Pre(PRN)%HisRange=eoshift(STA%STA(k)%Pre(PRN)%HisRange,shift=1,boundary=99.d0)
            STA%STA(k)%Pre(PRN)%HisPhase=eoshift(STA%STA(k)%Pre(PRN)%HisPhase,shift=1,boundary=99.d0)
             if ((P2/=0.d0) .and. (L2/=0.d0)) then
                STA%STA(k)%Pre(PRN)%HisRange(60)=P2
                STA%STA(k)%Pre(PRN)%HisPhase(60)=L2
                Range=P2
                Phase=L2
            end if
        else
            write(*,*) 'Smooth_Combine incorrect, please check: ', Smooth_Combine
            pause
            stop
        end if
        i=ceiling(Smooth_Time/Interval)
        if (i>59) i=59
        if ( STA%STA(k)%Pre(PRN)%n >i)  STA%STA(k)%Pre(PRN)%n=i
        call mean(STA%STA(k)%Pre(PRN)%HisRange(60-i:59),i,STA%STA(k)%Pre(PRN)%Range)
        call mean(STA%STA(k)%Pre(PRN)%HisPhase(60-i:59),i,STA%STA(k)%Pre(PRN)%PrePhase)
        if ((Range/=0.d0) .and. (Phase/=0.d0)) then   ! Smooth
            STA%STA(k)%Pre(PRN)%n=STA%STA(k)%Pre(PRN)%n+1.d0 ! dsind(Ele)
            Range=Range*1.d0/STA%STA(k)%Pre(PRN)%n+(1.d0-1.d0/STA%STA(k)%Pre(PRN)%n)*(STA%STA(k)%Pre(PRN)%Range+(Phase-STA%STA(k)%Pre(PRN)%PrePhase) )  
        end if
        if (index(Smooth_Combine,"P1")/=0) then
            P1=Range
        elseif ( (index(Smooth_Combine,"P2")/=0)  .and. (freq_comb=='L2L3')) then
            P1=Range
        else
            P2=Range
        end if
        ! Due to the ionophere opposite in Range and Phase
!        if (STA%STA(k)%Pre(PRN)%n>3.d0) STA%STA(k)%Pre(PRN)%n=3.d0
    else 
        write(*,*) '***ERROR: can not recognise smooth combination ', Smooth_Combine
        pause
        stop
    end if

    STA%STA(k)%Pre(PRN)%Range=Range   ! 平滑后的伪距
    STA%STA(k)%Pre(PRN)%PrePhase=Phase   ! 前一个历元的相位
end subroutine

! ============= mean ==========
!  Get the mean value of A.
!  If element in A is zero, skip it.
! INPUTS:
!     A        array of A
!     N        size of A
! OUTPUT:
!     MA     mean value of A
! WRITTEN BY: Yize Zhang, zhyize@163.com
! =============================
subroutine mean(A, N, MA)
implicit none
    real(8) :: A(N), MA
    integer :: N
    integer :: i, k
    
    k=0
    MA=0.d0
    do i=1,N
        if (dabs(A(i)-99.d0)>0.001d0) then
            k=k+1
            MA=MA+1.d0/k*(A(i)-MA)
        end if 
    end do

end subroutine