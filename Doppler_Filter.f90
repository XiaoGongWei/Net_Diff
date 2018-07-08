! ==================  Doppler_Filter  ================
!
! PURPOSE:
!            Smooth the pseudo-range combined with phase observations
!            by using Doppler flitering method.
!
! INPUTS:
!          Sow                Sow
!          PRN               PRN of the satellite
!          P1,P2             Pseudorange observation in L1 & L2 band, in meters
!          D1,D2              Phase observation in L1 & L2 band, in meters
!          k                     station order
! 
! OUTPUT:
!          P1, P2            Smoothed pseudo-range
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ===============  End of Header  ==================

subroutine Doppler_Filter(Sow, PRN, P1, P2, Range, DP1, DP2, k)
use MOD_FileID
use MOD_Var
use MOD_STA
implicit none
    ! Intent
    real(8) :: Sow
    integer     :: PRN, k
    real(8) :: P1, P2,Range, DP1, DP2
    real(8) :: DP, DP_Previous, t1, t2

    if (DP1/=0.d0 .and. DP2/=0.d0) then
        DP=(DP1*c/f1+DP2*c/f2)/2.d0
    elseif (DP1/=0.d0) then
        DP=DP1*c/f1
    elseif(DP1/=0.d0) then
        DP=DP2*c/f2
    else
        return
    end if
    DP_Previous=STA%STA(k)%Pre(PRN)%Dp
    STA%STA(k)%Pre(PRN)%Dp=Dp
    if (DP_Previous/=0.d0) then
        Dp=(DP_Previous+DP)/2.d0   ! Average of previous epoch and current epoch
    end if
    if (STA%STA(k)%Pre(PRN)%Sow1<=0.d0 .or. STA%STA(k)%Pre(PRN)%Sow2<=0.d0 ) then
        STA%STA(k)%Pre(PRN)%Sow1=Sow   ! First epoch
        STA%STA(k)%Pre(PRN)%Sow2=Sow
    end if
    t1=Sow-STA%STA(k)%Pre(PRN)%Sow1
    t2=Sow-STA%STA(k)%Pre(PRN)%Sow2
    if (t1-10.d0<0.d0) then  ! Data gap no less than 10 second
        STA%STA(k)%Pre(PRN)%n1=STA%STA(k)%Pre(PRN)%n1+1.d0
    else
        STA%STA(k)%Pre(PRN)%n1=1.d0
    end if
    if (t2-10.d0<0.d0) then
        STA%STA(k)%Pre(PRN)%n2=STA%STA(k)%Pre(PRN)%n2+1.d0
    else
        STA%STA(k)%Pre(PRN)%n2=1.d0
    end if

    if (P1/=0.d0) then   ! How about ionosphere effect?
!        write(CSID,'(I3,F8.1,F10.2)',advance='no') PRN, STA%STA(k)%Pre(PRN)%n1,-P1+STA%STA(k)%Pre(PRN)%P1-DP*t1
        P1=P1/STA%STA(k)%Pre(PRN)%n1+(1.d0-1.d0/STA%STA(k)%Pre(PRN)%n1)*(STA%STA(k)%Pre(PRN)%P1-DP*t1)
!    elseif (STA%STA(k)%Pre(PRN)%P1/=0.d0) then
!        P1=STA%STA(k)%Pre(PRN)%P1-DP*t1
!    end if
        STA%STA(k)%Pre(PRN)%Sow1=Sow  ! update time and smooth pseudo-range
    end if
    if (P2/=0.d0) then
!        write(CSID,'(F10.2)',advance='no') -P2+STA%STA(k)%Pre(PRN)%P2-DP*t2
        P2=P2/STA%STA(k)%Pre(PRN)%n2+(1.d0-1.d0/STA%STA(k)%Pre(PRN)%n2)*(STA%STA(k)%Pre(PRN)%P2-DP*t2)
!    elseif (STA%STA(k)%Pre(PRN)%P2/=0.d0) then
!        P2=STA%STA(k)%Pre(PRN)%P2-DP*t2
!    end if
        STA%STA(k)%Pre(PRN)%Sow2=Sow
    end if

!    write(CSID,'(A)') ''
    STA%STA(k)%Pre(PRN)%P1=P1
    STA%STA(k)%Pre(PRN)%P2=P2
    if ( (P1/=0.d0) .and. (P2/=0.d0) ) then
        Range=(f1*f1*P1-f2*f2*P2)/(f1+f2)/(f1-f2) ! Ionospheric-free combination
    end if

    return

end subroutine